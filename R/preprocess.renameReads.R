#' @title renameReads
#'
#' @description Renames and copies raw fastq.gz read files to a new directory
#'   using a sample spreadsheet that maps old names to new names. Supports
#'   samples spread across multiple lanes by specifying lane label strings via
#'   lane.names. Output files are named using the convention
#'   NewName_L00N_READ1/2.fastq.gz.
#'
#' @param sample.spreadsheet path to a CSV file with columns Old_Name and
#'   New_Name specifying the renaming mapping.
#'
#' @param output.directory path to the directory where renamed files will be
#'   copied.
#'
#' @param skip.not.found logical; if TRUE samples in the spreadsheet whose
#'   files cannot be located are silently skipped rather than raising an error.
#'
#' @param lane.names character vector of lane label substrings (e.g. c("s2",
#'   "s3")) used to distinguish reads from different lanes when a sample has
#'   multiple lanes; NULL assumes a single lane per sample.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before copying files.
#'
#' @return invisibly; side effect is renamed fastq.gz files written to
#'   output.directory.
#'
#' @export

renameReads = function(sample.spreadsheet = NULL,
                       output.directory = NULL,
                       skip.not.found = FALSE,
                       lane.names = NULL,
                       overwrite = FALSE){

  # ### Example usage
  # sample.spreadsheet = "/Volumes/Main_Data/Sequence_Capture/HF01_Ranoidea-Cornufer-twolanes_July2017/rename.csv"
  # output.directory = "/Volumes/Main_Data/Sequence_Capture/HF01_Ranoidea-Cornufer-renamed_July2017"
  # read.directory = "/Volumes/Main_Data/Sequence_Capture/HF01_Ranoidea-Cornufer-twolanes_July2017"
  # overwrite = FALSE
  # lane.names = c("s2", "s3")

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }

    if (overwrite == FALSE){ stop("Overwrite == FALSE and directory exists.") }
  }#end else

  sample.data = read.csv(sample.spreadsheet)
  all.reads = list.files(read.directory)
  all.reads = all.reads[grep("fastq.gz$|fq.gz$", all.reads)]

  if (nrow(sample.data) == 0){ return("no samples remain to analyze.") }

  for (i in 1:nrow(sample.data)){

    sample.reads = all.reads[grep(as.character(sample.data$Old_Name[i]), all.reads)]

    if (length(sample.reads) == 0){ stop("here") }

    if (is.null(lane.names) != T){

      for (j in 1:length(lane.names)){

        lane.reads = sample.reads[grep(paste0("_", lane.names[j], "_"), sample.reads)]

        #Save the read files with the new names in the new directory
        system(paste0("cp ", read.directory, "/", lane.reads[1], " ",
                      paste0(output.directory, "/", sample.data$New_Name[i], "_L00", j, "_READ1.fastq.gz")))

        system(paste0("cp ", read.directory, "/", lane.reads[2], " ",
                      paste0(output.directory, "/", sample.data$New_Name[i], "_L00", j, "_READ2.fastq.gz")))

      }#end j loop

    } else{

      if (length(sample.reads) >= 3){ stop("There are multiple lanes. Please add in laberls to lane.names")}

      #Save the read files with the new names in the new directory
      system(paste0("cp ", read.directory, "/", sample.reads[1], " ",
                    paste0(output.directory, "/", sample.data$New_Name[i], "_L001_READ1.fastq.gz")))

      system(paste0("cp ", read.directory, "/", sample.reads[2], " ",
                    paste0(output.directory, "/", sample.data$New_Name[i], "_L001_READ2.fastq.gz")))

    }#end if



  }#end i loop

}#end function

