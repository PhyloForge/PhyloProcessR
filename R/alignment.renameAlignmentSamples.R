#' @title renameAlignmentSamples
#'
#' @description Renames sample sequences within phylip alignment files in a directory (and optionally sub-directories) using a two-column lookup table. Each alignment is read, any sample names found in the lookup table are replaced with the corresponding new name, and the file is overwritten in place.
#'
#' @param align.directory path to the directory containing alignment files to be renamed
#'
#' @param align.extension file extension used to identify alignment files (e.g. ".phy")
#'
#' @param rename.file path to a tab-delimited text file with a header row; must contain columns named "Old_Name" and "New_Name" specifying the current and replacement sample names
#'
#' @param recursive if TRUE, search align.directory recursively to rename alignments in sub-directories as well; if FALSE, only rename alignments directly in align.directory
#'
#' @return overwrites each alignment file in place with the renamed sequences; nothing is returned to R
#'
#' @export

###Rename alignments
renameAlignmentSamples = function(align.directory = NULL,
                                  align.extension = ".phy",
                                  rename.file = NULL,
                                  recursive = TRUE){

  #Debug
  # align.directory = tree.dir
  # align.extension = ".phy"
  # rename.file = rename.path
  # recursive = TRUE

  #Reads in rename data
  rename.table = read.table(rename.file, header = T)

  #Reads in the files
  file.names = list.files(align.directory,
                          full.names = T,
                          recursive = recursive,
                          pattern = paste0(align.extension, "$") )

  ## Replace old names with new ones
  for(i in 1:length(file.names)){

    align = Biostrings::readDNAMultipleAlignment(file = file.names[i], format = "phylip")
    align = Biostrings::DNAStringSet(align)

    #Fix the labels
    temp.labels = names(align)
    temp.rename = rename.table[rename.table[,1] %in% temp.labels,]
    if (nrow(temp.rename) == 0){ next }
    temp.labels[pmatch(temp.rename[,1], temp.labels)] = as.character(temp.rename[,2])
    names(align) = temp.labels

    #Saves the file finally
    write.temp = strsplit(as.character(align), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
    writePhylip(aligned.set, file = file.names[i], interleave = F)

  }#end i loop

}#end function



