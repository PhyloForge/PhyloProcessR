#' @title makeFileRename
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param read.directory path to a folder of raw reads in fastq format.
#'
#' @param output.dir the new directory to save the adaptor trimmed sequences
#'
#' @param rename.file "Sample" to run on a single sample or "Directory" to run on a directory of samples
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
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
    } else {
      stop("File exists and overwrite = FALSE.")
    }
  }#end else

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
