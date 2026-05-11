#' @title trimAlignmentColumns
#'
#' @description Removes columns from a multiple sequence alignment where the proportion of gap characters ("-"), and optionally N and "?" characters, meets or exceeds a specified threshold. Alignments with two or fewer sequences are returned unmodified.
#'
#' @param alignment a DNAStringSet (or compatible format) containing the aligned sequences
#'
#' @param min.gap.percent percentage threshold (0-100): columns where gap/N/? characters comprise at least this percentage of sequences are removed
#'
#' @param count.n if TRUE, count "N" and "?" characters as gaps when computing the gap percentage; if FALSE, only "-" and "?" are counted
#'
#' @return a DNAStringSet with high-gap columns removed
#'
#' @export

trimAlignmentColumns = function(alignment = NULL,
                                min.gap.percent = 50,
                                count.n = T){
  #Debug
  # alignment = sim.align
  # min.gap.percent = 55

  if (length(alignment) <= 2){ return(alignment) }

  #Convert alignment to easier to work with matrix
  temp.align = strsplit(as.character(alignment), "")
  mat.align = lapply(temp.align, tolower)
  m.align = as.matrix(ape::as.DNAbin(mat.align))

  del.col = c()
  for (k in 1:ncol(m.align)){
    #Tables gap characters and calculataes prpoption
    gaps = table(as.character(m.align[,k]))

    #Counts Ns as gaps
    if (count.n == T){ per.gaps = sum(gaps[names(gaps) %in% c("-", "n", "?")]/nrow(m.align)) }

    #Doesn't count N as gaps
    if (count.n == F){ per.gaps = sum(gaps[names(gaps) %in% c("-", "?")]/nrow(m.align)) }

    if (length(per.gaps) == 0){ next }

    #Records column when the gaps exceed this percentage
    if (per.gaps * 100 >= min.gap.percent){ del.col = append(del.col, k) }

  }#end k loop

  #Removes bad columns and converts
  if (length(del.col) != 0){ del.align = m.align[,-del.col] } else { del.align = m.align }
  out.align = alignmentConversion(input.alignment = del.align,
                                  end.format = "DNAStringSet")

  return(out.align)

}#end function
