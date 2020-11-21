#' @title replaceAlignmentCharacter
#'
#' @description Function for replacing a character with another in a multiple sequence alignment
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param char.find the character to find
#'
#' @param char.replace the character to replace "char.find" with
#'
#' @return returns DNAStringSet of alignment with replaced characters
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export


replaceAlignmentCharacter = function(alignment = NULL,
                                     char.find = "N",
                                     char.replace = "-"){

  #applies across DNASTringSet. Convert formats?
  t.align = lapply(alignment, function (x) gsub("N", "-", x))
  n.align = Biostrings::DNAStringSet(unlist(t.align))

  return(n.align)

}#end function


