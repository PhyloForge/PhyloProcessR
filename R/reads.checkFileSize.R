#' @title reads.checkFileSize
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param input.reads path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param decontamination.path directory of genomes contaminants to scan samples
#'
#' @param samtools.path system path to samtools in case it can't be found
#'
#' @param bwa.path system path to bwa in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
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

reads.checkFileSize = function(read.directory = NULL,
                               check.file = NULL,
                               output.name = "read-fileSize",
                               overwrite = FALSE) {

  #Debug
  # setwd("/Volumes/Main_Data/Raw_Data/HF11_All-Frogs_Sept2021")
  # read.directory = "/Volumes/Main_Data/Raw_Data/HF11_All-Frogs_Sept2021/raw_data"
  # output.name = "read-counts-verified"
  # check.file = "checkSize.csv"
  # overwrite = TRUE

  #Quick checks
  if (is.null(read.directory) == TRUE){ stop("Please provide input reads.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.name) == TRUE && overwrite == TRUE){
      system(paste0("rm ", output.name))
  }#end else

  #Read in sample data **** sample is run twice?!
  reads = list.files(read.directory, recursive = T, full.names = T)
  reads = reads[grep("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", reads)]

  check.size = read.csv(check.file, header = FALSE)

  #Creates the summary log
  summary.data =  data.frame(File_Name = as.character(),
                             Read_fileSize = as.numeric(),
                             Read_checkSize = as.numeric(),
                             checkSize_Pass = as.character())

  #Runs through each sample
  for (i in 1:length(reads)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.reads = reads[i]

    sample.nm = gsub(".*\\/", "", sample.reads)
    sample.check = check.size[grep(sample.nm, check.size[,2]),]

    measure.size = file.size(sample.reads)
    actual.size = sample.check[1,1]

    pass.fail = "PASS"
    if (actual.size != measure.size){ pass.fail = "FAIL" }

    temp.remove = data.frame(File_Name = reads[i],
                               Read_fileSize = measure.size,
                               Read_checkSize = actual.size,
                               checkSize_Pass = pass.fail)

    summary.data = rbind(summary.data, temp.remove)

    print(paste0(reads[i], " Completed fastq counting!"))

  }#end sample i loop

  write.csv(summary.data, file = paste0(output.name, ".csv"), row.names = FALSE)
  return(summary.data)
}

