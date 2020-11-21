#' @title renameAlignmentSamples
#'
#' @description Function for easily renaming the samples of a multiple sequence alignment within a directory, or a single alignment
#'
#' @param align.directory the directory of alignments to have their samples renamed
#'
#' @param align.extension the file extension of your alignment files to be named
#'
#' @param rename.file a tab delimited text file with two columns named Old_Name and New_Name for the current and new name of the target sample
#'
#' @param recursive TRUE to recursively rename within sub-directories or FALSE just for the main directory
#'
#' @return overwrites alignments in directory with renamed sample alignments
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

    align = Biostrings::readAAMultipleAlignment(file = file.names[i], format = "phylip")
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



