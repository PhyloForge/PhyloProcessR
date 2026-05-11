#' @title utility.checkFileSize
#'
#' @description Verifies the integrity of raw read files by comparing each
#'   file's actual size on disk against expected sizes recorded in a reference
#'   size file. Produces a PASS/FAIL report per file and writes the results to
#'   a CSV.
#'
#' @param read.directory path to the directory containing the fastq.gz read
#'   files to be checked.
#'
#' @param check.file.name path to a two-column whitespace-delimited text file
#'   where column 1 is the expected file size in bytes and column 2 is the file
#'   name. Typically produced by running \code{ls -l} or a similar command on
#'   the source server.
#'
#' @param output.name base name (without extension) for the output CSV
#'   containing PASS/FAIL results.
#'
#' @param overwrite logical; if TRUE an existing output CSV is deleted before
#'   writing.
#'
#' @return a data frame with columns File_Name, Read_fileSize (measured),
#'   Read_checkSize (expected), and checkSize_Pass (PASS or FAIL). Also written
#'   to output.name.csv.
#'
#' @export

utility.checkFileSize = function(read.directory = NULL,
                                check.file.name = NULL,
                                output.name = "check-fileSize",
                                overwrite = FALSE) {

  #Debug
  # read.directory = "/Users/chutter/Dropbox/Research/3_Sequence-Database/Raw-Reads/HF13_Ancient_DNA_2/Reads"
  # output.name = "read-fileSize"
  # check.file = "/Users/chutter/Dropbox/Research/3_Sequence-Database/Raw-Reads/HF13_Ancient_DNA_2/checkSize.txt"
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

  check.size = read.table(check.file, header = FALSE)

  #Creates the summary log
  summary.data =  data.frame(File_Name = as.character(),
                            Read_fileSize = as.numeric(),
                            Read_checkSize = as.numeric(),
                            checkSize_Pass = as.character())

  #Runs through each sample
  for (i in seq_along(reads)) {
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

