# Internal helper: lists all files recursively in a Dropbox folder via API v2,
# avoiding the rdrop2::drop_dir() LinearizeNestedList bug on empty entries.
.dropbox_list_files = function(path, token) {
  resp = httr::POST(
    url = "https://api.dropboxapi.com/2/files/list_folder",
    httr::config(token = token),
    httr::content_type_json(),
    body = jsonlite::toJSON(
      list(path = path, recursive = TRUE, limit = 2000L),
      auto_unbox = TRUE
    )
  )
  httr::stop_for_status(resp)
  result = httr::content(resp, as = "parsed")
  entries = result$entries

  while (isTRUE(result$has_more)) {
    resp = httr::POST(
      url = "https://api.dropboxapi.com/2/files/list_folder/continue",
      httr::config(token = token),
      httr::content_type_json(),
      body = jsonlite::toJSON(list(cursor = result$cursor), auto_unbox = TRUE)
    )
    httr::stop_for_status(resp)
    result = httr::content(resp, as = "parsed")
    entries = c(entries, result$entries)
  }

  paths = vapply(entries, function(e) {
    if (!is.null(e[[".tag"]]) && e[[".tag"]] == "file") e$path_display else NA_character_
  }, character(1))
  paths[!is.na(paths)]
}


# Internal helper: interactive authentication. The Dropbox application key and
# secret come from the DROPBOX_APP_KEY and DROPBOX_APP_SECRET environment
# variables. Register an application at https://www.dropbox.com/developers to
# get them. Supply dropbox.token instead to skip interactive authentication.
.drop_auth = function() {
  app.key <- Sys.getenv("DROPBOX_APP_KEY")
  app.secret <- Sys.getenv("DROPBOX_APP_SECRET")

  if (nchar(app.key) == 0 || nchar(app.secret) == 0) {
    stop("Interactive Dropbox authentication needs the DROPBOX_APP_KEY and ",
         "DROPBOX_APP_SECRET environment variables. Set them, or give a saved ",
         "token file with the dropbox.token argument.")
  }

  dropbox <- httr::oauth_endpoint(
    authorize = "https://www.dropbox.com/oauth2/authorize",
    access = "https://api.dropbox.com/oauth2/token"
  )
  dropbox_app <- httr::oauth_app("dropbox", app.key, app.secret)
  httr::oauth2.0_token(dropbox, dropbox_app, cache = TRUE)
}

# Internal helper: download file
.drop_download = function(path, local_path, token, overwrite = FALSE) {
  if (file.exists(local_path) && !overwrite) stop("File exists")
  resp <- httr::POST(
    url = "https://content.dropboxapi.com/2/files/download",
    httr::config(token = token),
    httr::add_headers(
      `Dropbox-API-Arg` = jsonlite::toJSON(list(path = path), auto_unbox = TRUE)
    ),
    httr::write_disk(local_path, overwrite = overwrite)
  )
  httr::stop_for_status(resp)
  return(TRUE)
}


# Internal helper: puts a read pair in first-mate, second-mate order. Sorting on
# its own is not enough when the mate label is not the last part of the name.
.orderReadPair = function(read.files = NULL) {

  base.names = basename(read.files)
  first.mate = grepl("_1\\.f|-1\\.f|_R1|-R1|_READ1|-READ1", base.names)
  second.mate = grepl("_2\\.f|-2\\.f|_R2|-R2|_READ2|-READ2", base.names)

  if (sum(first.mate) == 1 && sum(second.mate) == 1) {
    return(c(read.files[first.mate], read.files[second.mate]))
  }

  return(sort(read.files))
}#end .orderReadPair


#' @title dropboxDownload
#'
#' @description Downloads paired-end fastq.gz read files from a Dropbox account
#'   using the Dropbox API v2. A sample spreadsheet maps source file names to
#'   desired sample names and the function renames files to the standard
#'   convention (SampleName_L00N_READ1/2.fastq.gz) on download. A
#'   file_rename_dropbox.csv is written on completion for use in downstream
#'   functions.
#'
#' @param sample.spreadsheet path to a CSV file with at least two columns: File
#'   (the file name prefix in Dropbox to search for) and Sample (the desired
#'   output sample name). Multiple rows per sample are treated as separate
#'   sequencing lanes.
#'
#' @param dropbox.directory the Dropbox path (starting with /) to the directory
#'   to search for read files.
#'
#' @param dropbox.token path to an RDS file containing a saved Dropbox OAuth2
#'   token. If NULL, the function starts interactive authentication, which
#'   needs the DROPBOX_APP_KEY and DROPBOX_APP_SECRET environment variables.
#'
#' @param output.directory local path where downloaded files will be saved.
#'
#' @param skip.not.found logical; if TRUE samples whose files cannot be located
#'   in Dropbox are silently skipped. If FALSE an error is raised.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before downloading. FALSE resumes and downloads only the lanes
#'   that are not yet complete.
#'
#' @return invisibly returns the rename table; writes downloaded fastq.gz files
#'   to output.directory and a file_rename_dropbox.csv in the working
#'   directory.
#'
#' @export

dropboxDownload = function(sample.spreadsheet = NULL,
                          dropbox.directory = NULL,
                          dropbox.token = NULL,
                          output.directory = NULL,
                          skip.not.found = FALSE,
                          overwrite = FALSE){

  #Quick checks
  if (is.null(sample.spreadsheet) == TRUE){ stop("Please provide a sample spreadsheet.") }
  if (file.exists(sample.spreadsheet) == F){ stop("Sample spreadsheet not found.") }
  if (is.null(dropbox.directory) == TRUE){ stop("Please provide a dropbox directory.") }
  if (is.null(output.directory) == TRUE){ stop("Please provide an output directory.") }
  if (is.null(dropbox.token) == FALSE){
    if (file.exists(dropbox.token) == F){ stop("Dropbox token file not found.") }
  }

  #Sets up the output directory
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE) { .resetDirectory(output.directory) }
  } # end else

  token = if (!is.null(dropbox.token)) readRDS(dropbox.token) else .drop_auth()
  all.reads = .dropbox_list_files(dropbox.directory, token)

  all.reads = all.reads[grep("fastq.gz$|fq.gz$", all.reads)]
  all.names = basename(all.reads)

  sample.data = read.csv(sample.spreadsheet)
  if (nrow(sample.data) == 0){ return("no samples available to download.") }

  sample.names = unique(sample.data$Sample)
  new.sample.data = data.frame(File = as.character(), Sample = as.character())
  for (i in seq_along(sample.names)){

    temp.data = sample.data[sample.data$Sample %in% sample.names[i], ]

    for (j in 1:nrow(temp.data)) {

      out.name = .sanitizeName(temp.data$Sample[j])
      lane.tag = sprintf("L%03d", j)
      outread.1 = paste0(output.directory, "/", out.name, "_", lane.tag, "_READ1.fastq.gz")
      outread.2 = paste0(output.directory, "/", out.name, "_", lane.tag, "_READ2.fastq.gz")

      # Skips a lane only when both read files are present and hold data. An
      # earlier version checked READ1 alone, so a download that stopped between
      # the two mates was never finished.
      if (.laneComplete(c(outread.1, outread.2)) == TRUE) {
        temp.sample.data <- data.frame(File = paste0(out.name, "_", lane.tag), Sample = out.name)
        new.sample.data <- rbind(new.sample.data, temp.sample.data)
        next
      }

      sample.reads = .matchPrefix(all.reads, all.names, as.character(temp.data$File[j]))

      # Checks the Sample column in case already renamed
      if (length(sample.reads) == 0) {
        sample.reads = .matchPrefix(all.reads, all.names, as.character(temp.data$Sample[j]))
      }

      # Last resort: the file name holds the search string somewhere else
      if (length(sample.reads) == 0) {
        sample.reads = all.reads[grep(as.character(temp.data$File[j]), all.names, fixed = TRUE)]
      }

      # For checking if reads are present
      if (length(sample.reads) >= 3) {
        stop("Problem with reads, sample ", sample.names[i], " File column matches to more than 1 sample. Check to ensure sample spreadsheet has multiple entries for samples with more than 1 lane of data.")
      }

      # Skip not found or crash
      if (length(sample.reads) == 0) {
        if (skip.not.found == FALSE) {
          stop(paste0("Error: sample reads for ", temp.data$Sample[j], " not found!"))
        } else {
          next
        }
      } # end if

      if (length(sample.reads) == 1) {
        if (skip.not.found == FALSE) {
          stop(paste0("Error: only one read set found for ", temp.data$Sample[j], " found!"))
        } else {
          next
        }
      } # end if

      # Save the read files with the new names in the new directory
      sample.reads = .orderReadPair(sample.reads)

      .drop_download(token = token,
        path = sample.reads[1],
        local_path = outread.1,
        overwrite = TRUE
      )

      .drop_download(token = token,
        path = sample.reads[2],
        local_path = outread.2,
        overwrite = TRUE
      )

      temp.sample.data = data.frame(File = paste0(out.name, "_", lane.tag), Sample = out.name)
      new.sample.data = rbind(new.sample.data, temp.sample.data)

    }#end j loop

  }#end i loop

  # The CSV is quoted so that a sample name with a comma cannot break the table
  write.csv(new.sample.data,
    file = "file_rename_dropbox.csv",
    row.names = FALSE
  )

  return(invisible(new.sample.data))
}#end function
