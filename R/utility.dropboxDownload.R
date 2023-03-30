#' @title dropboxDownload
#'
#' @description Function for downloading data off your dropbox account
#'
#' @param sample.spreadsheet spreadsheet path or table of file names to download
#'
#' @param dropbox.directory the directory that contains these target files
#'
#' @param dropbox.token your dropbox token path, created the first time you run the function but needs interactivity for log-in. Can be done once on your computer and the token can be moved around.
#'
#' @param output.directory the directory to save the files you are downloading
#'
#' @param skip.not.found the directory to save the files you are downloading
#'
#' @param overwrite whether to overwrite or not
#'
#' @return saves dropbox files to output directory
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export

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

  all.reads = rdrop2::drop_dir(
    path = dropbox.directory,
    recursive = TRUE
  )$path_display

  all.reads = all.reads[grep("fastq.gz$|fq.gz$", all.reads)]

  sample.data = read.csv(sample.spreadsheet)

  sample.names = unique(sample.data$Sample)
  new.sample.data = data.frame(File = as.character(), Sample = as.character())
  for (i in seq_along(sample.names)){

    temp.data = sample.data[sample.data$Sample %in% sample.names[i], ]
    
    for (j in 1:nrow(temp.data)) {
    
      sample.reads = all.reads[grep(as.character(temp.data$File[j]), all.reads)]

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
    file = "file_rename.csv",
    row.names = FALSE, quote = FALSE
  )

}#end function