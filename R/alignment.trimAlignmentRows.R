#' @title trimAlignmentRows
#'
#' @description Removes individual sequences (rows) from a multiple sequence alignment where the proportion of gap characters meets or exceeds a specified threshold. Alignments with two or fewer sequences are returned unmodified.
#'
#' @param alignment a DNAStringSet (or compatible format) containing the aligned sequences
#'
#' @param min.gap.percent percentage threshold (0-100): sequences where gap/N/? characters comprise at least this percentage of their length are removed
#'
#' @param count.n if TRUE, count "N" and "?" characters as gaps when computing the gap percentage per sequence; if FALSE, only "-" and "?" are counted
#'
#' @return a DNAStringSet with high-gap sequences removed
#'
#' @export

trimAlignmentRows = function(alignment = NULL,
                             min.gap.percent = 50,
                             count.n = T){
  #Debug
  # alignment = leg.align
  # min.gap.percent = 99
  # count.n = T

  if (length(alignment) <= 2){ return(alignment) }

  #Convert alignment to easier to work with matrix
  temp.align = strsplit(as.character(alignment), "")
  mat.align = lapply(temp.align, tolower)

  del.row = c()
  for (k in 1:length(mat.align)){
    #Tables gap characters and calculataes prpoption
    gaps = table(as.character(mat.align[[k]]))

    #Counts Ns as gaps
    if (count.n == T){ per.gaps = sum(gaps[names(gaps) %in% c("-", "n", "?")])/length(mat.align[[k]]) }

    #Doesn't count N as gaps
    if (count.n == F){ per.gaps = sum(gaps[names(gaps) %in% c("-", "?")])/length(mat.align[[k]]) }

    if (length(per.gaps) == 0){ next }

    #Records column when the gaps exceed this percentage
    if (per.gaps * 100 >= min.gap.percent){ del.row = append(del.row, k) }

  }#end k loop


  #Removes bad columns and converts
  if (length(del.row) != 0){ del.align = mat.align[-del.row] } else { del.align = mat.align }

  m.align = as.matrix(ape::as.DNAbin(del.align))

  out.align = PhyloCap::alignmentConversion(input.alignment = m.align,
                                            end.format = "DNAStringSet")

  return(out.align)

}#end function
