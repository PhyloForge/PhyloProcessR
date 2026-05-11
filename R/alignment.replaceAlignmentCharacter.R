#' @title replaceAlignmentCharacter
#'
#' @description Replaces all occurrences of a specified character with another character in every sequence of a DNAStringSet alignment. A common use is converting N characters to gap characters ("-") before trimming.
#'
#' @param alignment a DNAStringSet containing the aligned sequences to modify
#'
#' @param char.find the character to search for and replace in all sequences
#'
#' @param char.replace the character to substitute in place of char.find
#'
#' @return a DNAStringSet with all occurrences of char.find replaced by char.replace
#'
#' @export


replaceAlignmentCharacter = function(alignment = NULL,
                                     char.find = "N",
                                     char.replace = "-"){

  #applies across DNASTringSet. Convert formats?
  t.align = lapply(alignment, function (x) gsub("N", "-", x))
  n.align = Biostrings::DNAStringSet(unlist(t.align))

  return(n.align)

}#end function


