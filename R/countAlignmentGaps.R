#' @title countAlignmentGaps
#'
#' @description Function for counting the total number of gaps across the entire alignment
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @return returns a named vector of the number of gaps, total number of basepairs, and percent gaps in an entire alignment
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export

countAlignmentGaps = function(alignment = NULL){

  #alignment = non.align
  #Count gaps

  if (length(alignment) <= 2){
    gap.data = c(0, 100, 0)
    names(gap.data) = c("gaps", "alignment", "percent_gaps")
  }#end if

  len.temp = lapply(alignment, function (x) as.character(x))
  gap.temp = lapply(strsplit(as.character(len.temp), split = ""), function (x) length(which(x == "-")) )
  gap.count = unlist(gap.temp)
  total.temp = unlist(lapply(strsplit(as.character(len.temp), split = ""), function (x) length(x) ) )
  n.temp = lapply(strsplit(as.character(len.temp), split = ""), function (x) length(which(x == "N")) )
  q.temp = lapply(strsplit(as.character(len.temp), split = ""), function (x) length(which(x == "?")) )
  gap.count = gap.count + unlist(n.temp) + unlist(q.temp)
  gap.per = sum(gap.count)/sum(total.temp)

  if (is.nan(gap.per) == T){ gap.per = 0 }
  if (is.null(gap.per) == T){ gap.per = 0 }
  if (is.na(gap.per) == T){ gap.per = 0 }

  gap.data = c(sum(gap.count), sum(total.temp), gap.per * 100)
  names(gap.data) = c("gaps", "alignment", "percent_gaps")
  return(gap.data)

}
