#' @title makeFileRename
#'
#' @description Scans a directory of raw fastq files and generates a CSV
#'   rename table mapping file names to sample names. The resulting CSV can be
#'   supplied to organizeReads() or renameReads() to standardise file naming.
#'   Duplicate file entries are detected and flagged as errors.
#'
#' @param read.directory path to a directory containing raw fastq or fastq.gz
#'   files, either flat or in per-sample sub-folders.
#'
#' @param output.name file name (including .csv extension) for the rename table
#'   to be written.
#'
#' @param name.samples character; either "folder" to derive sample names from
#'   the immediate parent folder of each read file, or "file" to derive sample
#'   names by stripping read-pair suffixes (_R1, _R2, etc.) from the file name
#'   itself.
#'
#' @param subfolder logical; if TRUE and name.samples is "folder" or "file",
#'   the leading folder component of the path is stripped before extracting the
#'   file name, allowing one level of sub-directory nesting.
#'
#' @param overwrite logical; if TRUE an existing output.name CSV is deleted
#'   before writing.
#'
#' @return invisibly; writes a two-column CSV with columns File and Sample to
#'   output.name.
#'
#' @export

makeFileRename = function(read.directory = NULL,
                          output.name = "file_rename.csv",
                          name.samples = c("folder", "file"),
                          subfolder = FALSE,
                          overwrite = FALSE) {

  #Debegging
  # setwd("/Volumes/Rodents/Murinae")
  # read.directory = "/Volumes/Rodents/Murinae/Exome"
  # output.name = "file_rename.csv"
  # overwrite = FALSE
  # name.samples = "folder"
  # subfolder = TRUE

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(read.directory) == TRUE){ stop("Please provide a directory of raw reads.") }
  if (file.exists(read.directory) == F){ stop("Input reads not found.") }
  if (is.null(output.name) == TRUE){ stop("Please provide an output name.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.name) == T) {
    if (overwrite == TRUE){
      system(paste0("rm ", output.name))
    }#end else
  }

  reads = list.files(read.directory, recursive = T, full.names = F, )
  reads = reads[grep("fastq.gz$|fq.gz$|fastq$|fq$", reads)]

  #Folder is where each sample reads are in a folder
  if (name.samples == "folder"){

    if (subfolder == TRUE){
      file.names = gsub("^[^\\/]*\\/", "", reads)
    } else {
      file.names = gsub(".*\\/", "", reads)
    }#end else
    #obtains sample names
    sample.names = gsub("\\/.*", "", reads)
  }

  #Folder is where each sample reads are in a folder
  if (name.samples == "file"){
    #gathers names
    if (subfolder == TRUE){
      file.names = gsub("^[^\\/]*\\/", "", reads)
    } else {
      file.names = gsub(".*\\/", "", reads)
    }#end else
    sample.names = gsub(".*\\/", "", reads)
    sample.names = gsub("_R1.*", "", sample.names)
    sample.names = gsub("_R2.*", "", sample.names)
    sample.names = gsub("_READ1.*", "", sample.names)
    sample.names = gsub("_READ2.*", "", sample.names)
    sample.names = gsub("_1_.*", "", sample.names)
    sample.names = gsub("_2_.*", "", sample.names)
  }

  #Do some checks for duplicates and other weirdness
  file.rename = data.frame(File = file.names, Sample = sample.names)
  file.dupe = file.rename[duplicated(file.rename$File) == T,]

  if (nrow(file.dupe) != 0){
    print(file.dupe$File)
    stop("Error: The list of files printed are duplicated across different samples and/or lanes.")
  }

  #Gets file.names
  file.names = gsub("_R1.*", "", file.rename$File)
  file.names = gsub("_R2.*", "", file.names)
  file.names = gsub("_READ1.*", "", file.names)
  file.names = gsub("_READ2.*", "", file.names)
  file.names = gsub("_1_.*", "", file.names)
  file.names = gsub("_2_.*", "", file.names)
  file.names = gsub("_1.fastq.*", "", file.names)
  file.names = gsub("_2.fastq.*", "", file.names)
  file.rename$File = file.names

  final.file = file.rename[duplicated(file.rename$File) == F,]

  write.csv(final.file, file = output.name, row.names = F)

}#end function
