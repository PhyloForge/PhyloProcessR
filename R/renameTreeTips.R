#' @title renameTreeTips
#'
#' @description Function for easily renaming the tips of all phylogenetic trees within a directory (or a single tree)
#'
#' @param tree.directory the directory of trees to have their samples renamed
#'
#' @param tree.extension the file extension of your tree files to be named
#'
#' @param rename.file a tab delimited text file with two columns named Old_Name and New_Name for the current and new name of the target sample
#'
#' @param recursive TRUE to recursively rename within sub-directories or FALSE just for the main directory
#'
#' @return overwrites trees in directory with renamed sample trees
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

#function
renameTreeTips = function(tree.directory = NULL,
                          tree.extension = ".tre",
                          rename.file = NULL,
                          recursive = TRUE){

  #Debug
  # tree.directory = tree.dir
  # tree.extension = ".tre"
  # rename.file = rename.path
  # recursive = TRUE

  #Reads in rename data
  rename.table = read.table(rename.path, header = T)

  #Reads in the files
  file.names = list.files(tree.directory,
                          full.names = T,
                          recursive = recursive,
                          pattern = paste0(tree.extension, "$") )

  ## Replace old names with new ones
  for(i in 1:length(file.names)){

    temp.tree = ape::read.tree(file.names[i])

    if (class(temp.tree) == "multiPhylo"){
      #Goes through each tree
      for (j in 1:length(temp.tree)){
        temp.labels = temp.tree[[j]]$tip.label
        temp.rename = rename.table[rename.table[,1] %in% temp.labels,]
        temp.labels[pmatch(temp.rename[,1], temp.labels)] = as.character(temp.rename[,2])
        temp.tree[[j]]$tip.label = temp.labels
      }#end j loop
    }#end multiphylo

    #Single Tree
    if (class(temp.tree) == "phylo"){
      #Goes through each tree
      temp.labels = temp.tree$tip.label
      temp.rename = rename.table[rename.table[,1] %in% temp.labels,]
      temp.labels[pmatch(temp.rename[,1], temp.labels)] = as.character(temp.rename[,2])
      temp.tree$tip.label = temp.labels
    }#end multiphylo

    #Write the tree to file
    ape::write.tree(temp.tree, file = file.names[i])

  }#end i loop

}#end function

