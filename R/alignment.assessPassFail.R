#' @title alignmentAssess
#'
#' @description Function for assessing whether an entire alignment meets certain criteria
#'
#' @param alignment alignment in DNAStringSet Format
#'
#' @param max.alignment.gap.percent maximum percentage of gaps (including N and ?) permitted
#' across the whole alignment for it to pass. Alignments at or above this threshold fail.
#'
#' @param min.taxa.alignment minimum number of sequences required for the alignment to pass.
#' Alignments with this many or fewer sequences fail.
#'
#' @param min.alignment.length minimum alignment length (width in bp) required to pass.
#' Alignments at or below this length fail.
#'
#' @return a logical value, TRUE or FALSE, if your alignment passes (TRUE) or fails (FALSE) the filters
#'
#' @examples
#' \dontrun{
#' pass <- alignmentAssess(alignment = my.alignment,
#'                         max.alignment.gap.percent = 50,
#'                         min.taxa.alignment = 4,
#'                         min.alignment.length = 100)
#' }
#'
#' @export

#Assesses the alignment returning TRUE for pass and FALSE for fail
alignmentAssess = function(alignment = NULL,
                           max.alignment.gap.percent = 0,
                           min.taxa.alignment = 0,
                           min.alignment.length = 0){
  #Count gaps
  #alignment = non.align
  gap.data = countAlignmentGaps(alignment)
  #if (length(gap.data) != 3){ return(FALSE) }

  #Checks for enough data
  if (length(gap.data) <= 2){ return(FALSE) }

  #Checks for min percentage
  if (gap.data[3] >= max.alignment.gap.percent) {
    return(FALSE)
  }

  #Make a summary table thing
  if (length(alignment) <= min.taxa.alignment){
    return(FALSE)
  }

  if (max(Biostrings::width(alignment)) <= min.alignment.length){
    return(FALSE)
  }

  return(TRUE)

}#end function

