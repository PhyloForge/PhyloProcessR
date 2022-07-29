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

