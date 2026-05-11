#' @title organizeReads
#'
#' @description Copies raw fastq files from a source directory into a
#'   standardised per-sample directory structure, renaming files to a
#'   consistent convention (SampleName_L00N_READ1/2.fastq.gz). Sample-to-file
#'   mapping is provided by a two-column CSV (File, Sample) such as that
#'   produced by makeFileRename(). Each row in the CSV represents one
#'   paired-end lane for a sample; multiple rows for the same sample generate
#'   multiple lanes.
#'
#' @param read.directory path to the directory containing the raw fastq or
#'   fastq.gz files to be organised.
#'
#' @param output.directory path to the directory where the reorganised reads
#'   will be saved. A sub-directory is created for each unique sample name.
#'
#' @param rename.file path to a CSV file with at least two columns: File
#'   (partial or full file name used to locate the reads) and Sample (the
#'   desired output sample name).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated. Completed samples are skipped when overwrite is FALSE.
#'
#' @return invisibly; side effect is a populated output.directory with one
#'   sub-directory per sample containing renamed fastq.gz read pairs.
#'
#' @export

organizeReads = function(read.directory = NULL,
                        output.directory = "organized-reads",
                        rename.file = NULL,
                        overwrite = FALSE) {

  #Debegging
  # read.directory = "/Volumes/Backup_Hub/Raw_Data/HF14_Ultimate_FrogCap/Modern"
  # output.directory = "/Volumes/LaCie/Microhylidae/organized-reads"
  # rename.file = "/Users/chutter/Dropbox/Research/9_Paperwork/Arbor_Biosciences/HF-13/file_rename.csv"
  # overwrite = FALSE

  # #Debegging
  # read.directory = "/Volumes/Extreme_SSD/Unusual_Hosts/raw_reads"
  # output.directory = "/Volumes/Extreme_SSD/Unusual_Hosts/organized-reads"
  # rename.file = "/Volumes/Extreme_SSD/Unusual_Hosts/sample_names.csv"
  # overwrite = FALSE

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(read.directory) == TRUE){ stop("Please provide a directory of raw reads.") }
  if (file.exists(read.directory) == F){ stop("Input reads not found.") }
  if (is.null(rename.file) == TRUE){ stop("Please provide a table of file to sample name conversions.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Read in sample data and finds reads
  reads = list.files(read.directory, recursive = T, full.names = T)
  reads = reads[grep("fastq.gz$|fastq$|fq.gz$|fq$", reads)]

  sample.data = read.csv(rename.file)
  if (nrow(sample.data) == 0){ return("no samples available to organize.") }

  #Skips samples already finished
  if (overwrite == FALSE){
    done.names = list.files(output.directory)
    sample.data = sample.data[!sample.data$Sample %in% done.names,]
  } else { sample.data = sample.data }

  if (nrow(sample.data) == 0){ return("no samples available to organize.") }

  sample.names = unique(sample.data$Sample)

  for (i in seq_along(sample.names)){

    temp.data = sample.data[sample.data$Sample %in% sample.names[i], ]

    for (j in 1:nrow(temp.data)) {
      #################################################
      ### Part A: prepare for loading and checks
      #################################################
      # Finds all files for this given sample and turns into a cat string
      sample.reads = reads[grep(pattern = paste0(temp.data$File[j], "_"), x = reads)]
      # Checks the Sample column in case already renamed
      if (length(sample.reads) == 0) {
        sample.reads = reads[grep(pattern = paste0(temp.data$Sample[j], "_"), x = reads)]
      }
      if (length(sample.reads) == 0) {
        sample.reads = reads[grep(pattern = temp.data$File[j], x = reads)]
      }
      if (length(sample.reads) == 0) {
        sample.reads = reads[grep(pattern = temp.data$Sample[j], x = reads)]
      }

      # Returns an error if reads are not found
      if (length(sample.reads) == 0) {
        stop(paste0(
          temp.data$Sample[j], " does not have any reads present for files ",
          temp.data$File[j], " from the input spreadsheet."
        ))
      } # end if statement

      if (length(sample.reads) == 1) {
        stop(paste0(temp.data$Sample[j], " only one read file found. The other is missing."))
      }

      if (length(sample.reads) >= 3) {
        stop(paste0(temp.data$Sample[j], " had more than 2 read files associated, all file names must be unique."))
      }

      #################################################
      ### Part B: Create directories and move files
      #################################################
      # Create sample directory
      out.path = paste0(output.directory, "/", temp.data$Sample[j])
      if (file.exists(out.path) == FALSE) {
        dir.create(out.path)
      }

      # sets up output reads
      outread.1 = paste0(out.path, "/", temp.data$Sample[j], "_L00", j, "_READ1.fastq.gz")
      outread.2 = paste0(out.path, "/", temp.data$Sample[j], "_L00", j, "_READ2.fastq.gz")

      system(paste0("cp ", sample.reads[1], " ", outread.1))
      system(paste0("cp ", sample.reads[2], " ", outread.2))
    } # end j loop

  }#end i loop

}#end function
