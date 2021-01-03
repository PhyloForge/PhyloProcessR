#' @title convertAmbiguousConsensus
#'
#' @description Function for trimming out alignment columns with too many gaps
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param convert.r converts ambiguous base pair R to A or G
#'
#' @param convert.y converts ambiguous base pair Y to C or T
#'
#' @param convert.s converts ambiguous base pair S to G or C
#'
#' @param convert.w converts ambiguous base pair W to A or T
#'
#' @param convert.k converts ambiguous base pair K to G or T
#'
#' @param convert.m converts ambiguous base pair M to A or C
#'
#' @param convert.b converts ambiguous base pair B to C, G, or T
#'
#' @param convert.d converts ambiguous base pair D to A, G, or T
#'
#' @param convert.h converts ambiguous base pair H to A, C, or T
#'
#' @param convert.v converts ambiguous base pair V to A, C, or G
#'
#' @return returns DNAStringSet of column trimmed alignment
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export

convertAmbiguousConsensus = function(alignment = NULL,
                                     convert.r = c("A", "G"),
                                     convert.y = c("C", "T"),
                                     convert.s = c("G", "C"),
                                     convert.w = c("A", "T"),
                                     convert.k = c("G", "T"),
                                     convert.m = c("A", "C"),
                                     convert.b = c("C", "G", "T"),
                                     convert.d = c("A", "G", "T"),
                                     convert.h = c("A", "C", "T"),
                                     convert.v = c("A", "C", "G")) {
  #To do
  # Add randomization
  # Add select best
  #Add column?

  #Debug
 # alignment = non.align
  convert.r = "A"
  convert.y = "T"
  convert.s = "G"
  convert.w = "A"
  convert.k = "T"
  convert.m = "A"
  convert.b = "T"
  convert.d = "T"
  convert.h = "T"
  convert.v = "A"

  #applies across DNASTringSet. Convert formats?
  edit.align = lapply(alignment, function (x) gsub("R|r", convert.r, x))
  edit.align = lapply(edit.align, function (x) gsub("Y|y", convert.y, x))
  edit.align = lapply(edit.align, function (x) gsub("S|s", convert.s, x))
  edit.align = lapply(edit.align, function (x) gsub("W|w", convert.w, x))
  edit.align = lapply(edit.align, function (x) gsub("K|k", convert.k, x))
  edit.align = lapply(edit.align, function (x) gsub("M|m", convert.m, x))
  edit.align = lapply(edit.align, function (x) gsub("B|b", convert.b, x))
  edit.align = lapply(edit.align, function (x) gsub("D|d", convert.d, x))
  edit.align = lapply(edit.align, function (x) gsub("H|h", convert.h, x))
  edit.align = lapply(edit.align, function (x) gsub("V|v", convert.v, x))

  n.align = Biostrings::DNAStringSet(unlist(edit.align))

  return(n.align)

}#end function


