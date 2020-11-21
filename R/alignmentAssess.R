#' @title alignmentAssess
#'
#' @description Function for assessing whether an entire alignment meets certain criteria
#'
#' @param alignment alignment in DNAStringSet Format
#'
#' @param min.gap.percent minimum percentage of gaps allowed for an alignment to pass
#'
#' @param min.taxa.count minimum number of taxa allowed for an alignment to pass
#'
#' @param min.align.length minimum alignment length (width) allowed for an alignment to pass
#'
#' @return a logical value, TRUE or FALSE, if your alignment passes (TRUE) or fails (FALSE) the filters
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

#Assesses the alignment returning TRUE for pass and FALSE for fail
alignmentAssess = function(alignment = NULL,
                           min.gap.percent = 0,
                           min.taxa.count = 0,
                           min.align.length = 0){
  #Count gaps
  gap.data = countAlignmentGaps(alignment)

  #RECORD DATA
  if (gap.data[3] >= min.gap.percent) {
    return(FALSE)
  }

  #Make a summary table thing
  if (length(alignment) <= min.taxa.count){
    return(FALSE)
  }

  if (max(Biostrings::width(alignment)) <= min.align.length){
    return(FALSE)
  }

  return(TRUE)

}#end function

