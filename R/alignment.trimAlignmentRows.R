#' @title trimAlignmentRows
#'
#' @description Function for trimming out alignment rows (or samples) with too many gaps
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param min.gap.percent minimum threshold gap percentage allowed to trim row/sample
#'
#' @param count.n TRUE to count N as a gap
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
