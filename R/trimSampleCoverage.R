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
                              min.coverage.bp = 60,
                              relative.width = c("alignment", "sample")) {

  #Debug
  # alignment = non.align
  # min.coverage.bp = 300
  # min.coverage.percent = 50
  # relative.width = "sample"

  if (length(relative.width) == 2){ stop("please choose a relative width.") }

  if (length(alignment) <= 3){ return(alignment) }

  if (min.coverage.bp >= Biostrings::width(alignment)[1]){ return(alignment) }

  #Remove gap only alignments
  c.align = strsplit(as.character(alignment), "")
  gap.align = lapply(c.align, function(x) gsub("N|n", "-", x) )
  base.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  base.rem = base.count[base.count <= as.numeric(min.coverage.bp)]

  if (relative.width == "alignment"){
    r.width = Biostrings::width(alignment)[1]
  }

  if (relative.width == "sample"){
    r.width = max(base.count)
  }

  base.per = base.count/r.width
  base.rem = append(base.rem, base.per[base.per * 100 <= as.numeric(min.coverage.percent)])

  #Removes the bad stuff
  align.out = alignment[!names(alignment) %in% unique(names(base.rem))]
  return(align.out)

}#end function
