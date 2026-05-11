#' @title trimSampleCoverage
#'
#' @description Removes samples from an alignment whose sequence coverage falls below minimum thresholds. Coverage can be measured relative to the total alignment length ("alignment") or relative to the longest sample in the alignment ("sample"). Samples are removed if they have fewer than min.coverage.bp non-gap base pairs or if their non-gap base pair count is less than min.coverage.percent of the reference width. Alignments with three or fewer sequences or where min.coverage.bp exceeds the alignment length are returned unmodified.
#'
#' @param alignment a DNAStringSet containing the aligned sequences to filter
#'
#' @param min.coverage.percent minimum percentage of the reference width (alignment or longest sample) that each sample must cover in non-gap base pairs to be retained
#'
#' @param min.coverage.bp minimum absolute number of non-gap base pairs a sample must have to be retained
#'
#' @param relative.width whether the percentage threshold is calculated relative to "alignment" (total alignment width) or "sample" (the length of the longest sample in the alignment)
#'
#' @return a DNAStringSet with low-coverage samples removed
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
