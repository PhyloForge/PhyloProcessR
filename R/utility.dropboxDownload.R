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

  file.entries = entries[vapply(entries, function(entry) {
    is.null(entry[[".tag"]]) == FALSE && entry[[".tag"]] == "file"
  }, logical(1))]
  data.frame(
    Path = vapply(file.entries, function(entry) entry$path_display, character(1)),
    Size = vapply(file.entries, function(entry) {
      if (is.null(entry$size)) NA_real_ else as.numeric(entry$size)
    }, numeric(1)),
    ContentHash = vapply(file.entries, function(entry) {
      if (is.null(entry$content_hash)) NA_character_ else as.character(entry$content_hash)
    }, character(1)),
    stringsAsFactors = FALSE
  )
}


# Internal helper: interactive authentication. The Dropbox application key and
# secret come from the DROPBOX_APP_KEY and DROPBOX_APP_SECRET environment
# variables. Register an application at https://www.dropbox.com/developers to
# get them. Supply dropbox.token instead to skip interactive authentication.
.drop_auth = function() {
  app.key = Sys.getenv("DROPBOX_APP_KEY")
  app.secret = Sys.getenv("DROPBOX_APP_SECRET")

  if (nchar(app.key) == 0 || nchar(app.secret) == 0) {
    stop("Interactive Dropbox authentication needs the DROPBOX_APP_KEY and ",
         "DROPBOX_APP_SECRET environment variables. Set them, or give a saved ",
         "token file with the dropbox.token argument.")
  }

  dropbox = httr::oauth_endpoint(
    authorize = "https://www.dropbox.com/oauth2/authorize",
    access = "https://api.dropbox.com/oauth2/token"
  )
  dropbox.app = httr::oauth_app("dropbox", app.key, app.secret)
  httr::oauth2.0_token(dropbox, dropbox.app, cache = TRUE)
}

# Internal helper: download file
.drop_download = function(path, local.path, token, overwrite = FALSE) {
  if (file.exists(local.path) && !overwrite) stop("File exists")
  resp = httr::POST(
    url = "https://content.dropboxapi.com/2/files/download",
    httr::config(token = token),
    httr::add_headers(
      `Dropbox-API-Arg` = jsonlite::toJSON(list(path = path), auto_unbox = TRUE)
    ),
    httr::write_disk(local.path, overwrite = overwrite)
  )
  httr::stop_for_status(resp)
  return(TRUE)
}


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
  if (length(skip.not.found) != 1 || is.logical(skip.not.found) == FALSE || is.na(skip.not.found)) {
    stop("skip.not.found must be TRUE or FALSE.")
  }
  if (length(overwrite) != 1 || is.logical(overwrite) == FALSE || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  if (is.null(dropbox.token) == FALSE){
    if (file.exists(dropbox.token) == F){ stop("Dropbox token file not found.") }
  }

  sample.data = read.csv(sample.spreadsheet, stringsAsFactors = FALSE)
  sample.data = .validateRenameTable(sample.data, sanitize.samples = TRUE)
  if (nrow(sample.data) == 0){ return("no samples available to download.") }
  .checkFileOutsideOutput(sample.spreadsheet, output.directory)

  #Sets up the output directory
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE) { .resetDirectory(output.directory) }
  } # end else

  if (dir.exists("logs/sample_logs") == FALSE) { dir.create("logs/sample_logs", recursive = TRUE) }

  token = if (!is.null(dropbox.token)) readRDS(dropbox.token) else .drop_auth()
  dropbox.files = .dropbox_list_files(dropbox.directory, token)
  dropbox.files = dropbox.files[grep("fastq.gz$|fq.gz$", dropbox.files$Path), , drop = FALSE]
  all.reads = dropbox.files$Path
  all.names = basename(all.reads)

  sample.names = unique(sample.data$Sample)
  new.sample.data = data.frame(File = as.character(), Sample = as.character())
  for (i in seq_along(sample.names)){

    temp.data = sample.data[sample.data$Sample %in% sample.names[i], ]

    for (j in 1:nrow(temp.data)) {

      out.name = temp.data$Sample[j]
      lane.tag = sprintf("L%03d", j)
      outread.1 = paste0(output.directory, "/", out.name, "_", lane.tag, "_READ1.fastq.gz")
      outread.2 = paste0(output.directory, "/", out.name, "_", lane.tag, "_READ2.fastq.gz")

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
      # Skip not found or crash
      if (length(sample.reads) == 0) {
        if (skip.not.found == FALSE) {
          stop(paste0("Error: sample reads for ", temp.data$Sample[j], " not found!"))
        } else {
          next
        }
      } # end if

      if (length(sample.reads) != 2) {
        if (skip.not.found == FALSE) {
          stop("Expected exactly two read files for ", temp.data$Sample[j],
               ", but found ", length(sample.reads), ".")
        } else {
          next
        }
      } # end if

      # Save the read files with the new names in the new directory
      sample.reads = .orderReadPair(sample.reads)
      source.rows = match(sample.reads, dropbox.files$Path)
      metadata.file = file.path("logs/sample_logs", out.name,
                                paste0(out.name, "_", lane.tag, "_dropbox-metadata.csv"))
      metadata = c("source.1" = sample.reads[1], "source.2" = sample.reads[2],
                   "source.1.size" = dropbox.files$Size[source.rows[1]],
                   "source.2.size" = dropbox.files$Size[source.rows[2]],
                   "source.1.hash" = dropbox.files$ContentHash[source.rows[1]],
                   "source.2.hash" = dropbox.files$ContentHash[source.rows[2]],
                   "parameter.file" = temp.data$File[j],
                   "parameter.sample" = out.name)

      if (.laneComplete(c(outread.1, outread.2), metadata.file = metadata.file,
                        metadata = metadata) == TRUE) {
        temp.sample.data = data.frame(File = paste0(out.name, "_", lane.tag), Sample = out.name)
        new.sample.data = rbind(new.sample.data, temp.sample.data)
        next
      }
      if (overwrite == FALSE && .metadataConflicts(metadata.file, metadata) == TRUE) {
        stop(out.name, " ", lane.tag, " was downloaded from different source files. ",
             "Use overwrite = TRUE to replace it.")
      }

      temp.reads = c(tempfile(pattern = paste0(basename(outread.1), "-"),
                              tmpdir = output.directory, fileext = ".fastq.gz"),
                     tempfile(pattern = paste0(basename(outread.2), "-"),
                              tmpdir = output.directory, fileext = ".fastq.gz"))
      on.exit(unlink(temp.reads), add = TRUE)

      .drop_download(token = token,
        path = sample.reads[1],
        local.path = temp.reads[1],
        overwrite = TRUE
      )

      .drop_download(token = token,
        path = sample.reads[2],
        local.path = temp.reads[2],
        overwrite = TRUE
      )

      if (all(file.info(temp.reads)$size > 0) == FALSE) {
        stop("Dropbox returned an empty read file for ", out.name, " ", lane.tag, ".")
      }
      .publishFiles(temp.reads, c(outread.1, outread.2))
      .writeLaneMetadata(metadata, metadata.file)

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
