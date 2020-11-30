#' @title trimSampleCoverage
#'
#' @description Function for trimming out low coverage (i.e. high gap) samples from alignment
#'
#' @param alignment folder that contains an alignment in DNAStringSet format
#'
#' @param min.coverage.percent minimum percent sequence data present for each sample out of the total alignment length
#'
#' @param min.sample.bp minimum number of base pairs a sample needs to be kept
#'
#' @return returns sample coverage trimmed alignment. Alignment will be returned unmodified if the min.sample.bp is less than the alignment length.
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export


#Function trims for sample coverage
trimSampleCoverage = function(alignment = NULL,
                              min.coverage.percent = 50,
                              min.sample.bp = 60) {

  #Debug
  #alignment = edge.align
  #min.coverage.bp = 300
  #min.coverage.percent = 50

  if (length(alignment) <= 2){ return(alignment) }

  if (min.sample.bp >= width(alignment)[1]){ return(alignment) }

  #Remove gap only alignments
  c.align = strsplit(as.character(alignment), "")
  gap.align = lapply(c.align, function(x) gsub("N|n", "-", x) )
  gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  gap.rem = gap.count[gap.count <= as.numeric(min.sample.bp)]

  gap.count = gap.count/Biostrings::width(alignment)
  gap.rem = append(gap.rem, gap.count[gap.count * 100 <= as.numeric(min.coverage.percent)])

  #Removes the bad stuff
  align.out = alignment[!names(alignment) %in% unique(names(gap.rem))]
  return(align.out)

}#end function
