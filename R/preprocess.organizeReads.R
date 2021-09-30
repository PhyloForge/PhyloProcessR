#' @title organizeReads
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
#' @param resume TRUE to overwrite a folder of samples with output.dir
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

organizeReads = function(read.directory = NULL,
                         output.dir = "organized-reads",
                         rename.file = NULL,
                         overwrite = FALSE) {

  #Debegging
  # output.dir = "organized-reads"
  # rename.file = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/file_rename.csv"
  # overwrite = FALSE

  # setwd("/Volumes/Rodents/Murinae")
  # read.directory = "/Volumes/Rodents/Murinae/Exome"
  # dir.create("processed-reads")
  # output.dir = "processed-reads/organized-reads"
  # rename.file = "file_rename.csv"
  # overwrite = FALSE

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(read.directory) == TRUE){ stop("Please provide a directory of raw reads.") }
  if (file.exists(read.directory) == F){ stop("Input reads not found.") }
  if (is.null(rename.file) == TRUE){ stop("Please provide a table of file to sample name conversions.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.dir) == F){ dir.create(output.dir) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end else

  #Read in sample data and finds reads
  reads = list.files(read.directory, recursive = T, full.names = T)
  reads = reads[grep("fastq.gz$|fastq$|fq.gz$|fq$", reads)]

  sample.data = read.csv(rename.file)
  if (nrow(sample.data) == 0){ return("no samples available to organize.") }

  #Skips samples already finished
  if (overwrite == FALSE){
    done.names = list.files(output.dir)
    sample.data = sample.data[!sample.data$Sample %in% done.names,]
  } else { sample.data = sample.data }

  if (nrow(sample.data) == 0){ return("no samples available to organize.") }

  sample.count = 1
  for (i in 1:nrow(sample.data)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    #Finds all files for this given sample and turns into a cat string
    sample.reads = reads[grep(pattern = paste0(sample.data$File[i], "_"), x = reads)]
    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = paste0(sample.data$Sample[i], "_"), x = reads)] }
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.data$File[i], x = reads)] }
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.data$Sample[i], x = reads)] }

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(paste0(sample.data$Sample[i], " does not have any reads present for files ",
           sample.data$File[i], " from the input spreadsheet.") )
    } #end if statement

    if (length(sample.reads) == 1){
      stop(paste0(sample.data$Sample[i], " only one read file found. The other is missing."))
    }

    if (length(sample.reads) >= 3){
      stop(paste0(sample.data$Sample[i], " had more than 2 read files associated, all file names must be unique."))
    }

    #################################################
    ### Part B: Create directories and move files
    #################################################
    #Create sample directory
    out.path = paste0(output.dir, "/", sample.data$Sample[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path) }

    #Skips first iteration
    if (i != 1){
      if (sample.data$Sample[i] == sample.data$Sample[i-1]){
        sample.count = sample.count + 1
      } else {sample.count = 1 }
    }#end if

    #sets up output reads
    outread.1 = paste0(out.path, "/", sample.data$Sample[i], "_L00", sample.count, "_READ1.fastq.gz")
    outread.2 = paste0(out.path, "/", sample.data$Sample[i], "_L00", sample.count, "_READ2.fastq.gz")

    system(paste0("cp ", sample.reads[1], " ", outread.1))
    system(paste0("cp ", sample.reads[2], " ", outread.2))

  }#end i loop

}#end function
