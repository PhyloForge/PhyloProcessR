#' @title dropboxDownload
#'
#' @description Downloads paired-end fastq.gz read files from a Dropbox account
#'   using the rdrop2 package and the Dropbox API v2. A sample spreadsheet maps
#'   source file names to desired sample names and the function renames files to
#'   the standard convention (SampleName_L00N_READ1/2.fastq.gz) on download.
#'   A file_rename_dropbox.csv is written on completion for use in downstream
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
#'   token produced by rdrop2::drop_auth(). If NULL, the function attempts to
#'   retrieve a token from the rdrop2 cache.
#'
#' @param output.directory local path where downloaded files will be saved.
#'
#' @param skip.not.found logical; if TRUE samples whose files cannot be located
#'   in Dropbox are silently skipped. If FALSE an error is raised.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before downloading.
#'
#' @return invisibly; writes downloaded fastq.gz files to output.directory and
#'   a file_rename_dropbox.csv in the working directory.
#'
#' @export

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

dropboxDownload = function(sample.spreadsheet = NULL,
                          dropbox.directory = NULL,
                          dropbox.token = NULL,
                          output.directory = NULL,
                          skip.not.found = FALSE,
                          overwrite = FALSE){

  # # ### Example usage
  # setwd("/Volumes/LaCie/Bufonidae")
  # sample.spreadsheet <- "/Volumes/LaCie/Bufonidae/scripts/Bufonidae.csv"
  # output.directory <- "/Volumes/LaCie/Bufonidae/raw-reads"
  # dropbox.directory <- "/Research/3_Sequence-Database/Raw-Reads"
  # dropbox.token = "/Volumes/LaCie/dropbox-token.RDS"
  # rdrop2::drop_auth(rdstoken = dropbox.token)
  # dropbox.token <- "/Users/chutter/Dropbox/dropbox-token.RDS"
  # overwrite <- FALSE

  # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } # end else

  token = if (!is.null(dropbox.token)) readRDS(dropbox.token) else rdrop2:::get_dropbox_token()
  all.reads = .dropbox_list_files(dropbox.directory, token)

  all.reads = all.reads[grep("fastq.gz$|fq.gz$", all.reads)]

  sample.data = read.csv(sample.spreadsheet)

  sample.names = unique(sample.data$Sample)
  new.sample.data = data.frame(File = as.character(), Sample = as.character())
  for (i in seq_along(sample.names)){

    temp.data = sample.data[sample.data$Sample %in% sample.names[i], ]
    
    for (j in 1:nrow(temp.data)) {
    
      sample.reads = all.reads[grep(paste0(as.character(temp.data$File[j]), "_"), all.reads)]

      if (length(sample.reads) == 0) {
        sample.reads = all.reads[grep(paste0(as.character(temp.data$File[j]), "-"), all.reads)]
      }

      if (length(sample.reads) == 0) {
        sample.reads = all.reads[grep(as.character(temp.data$File[j]), all.reads, fixed = TRUE)]
      }

      if (file.exists(paste0(output.directory, "/", temp.data$Sample[j], "_L00", j, "_READ1.fastq.gz")) == TRUE) {
        temp.sample.data <- data.frame(File = paste0(temp.data$Sample[j], "_L00", j), Sample = temp.data$Sample[j])
        new.sample.data <- rbind(new.sample.data, temp.sample.data)
        next
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
      rdrop2::drop_download(
        path = sample.reads[1],
        local_path = paste0(output.directory, "/", temp.data$Sample[j], "_L00", j, "_READ1.fastq.gz"),
        overwrite = TRUE
      )

      rdrop2::drop_download(
        path = sample.reads[2],
        local_path = paste0(output.directory, "/", temp.data$Sample[j], "_L00", j, "_READ2.fastq.gz"),
        overwrite = TRUE
      )

      temp.sample.data = data.frame(File = paste0(temp.data$Sample[j], "_L00", j), Sample = temp.data$Sample[j])
      new.sample.data = rbind(new.sample.data, temp.sample.data)

    }#end j loop

  }#end i loop

  write.csv(new.sample.data,
    file = "file_rename_dropbox.csv",
    row.names = FALSE, quote = FALSE
  )

}#end function