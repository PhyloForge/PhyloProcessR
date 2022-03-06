#' @title trimAlignmentColumns
#'
#' @description Function for trimming out alignment columns with too many gaps
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param min.gap.percent minimum threshold gap percentage allowed to trim column
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

trimAlignmentColumns = function(alignment = NULL,
                                min.gap.percent = 100){
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
    per.gaps = sum(gaps[names(gaps) %in% c("-", "n", "?")]/nrow(m.align))

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
