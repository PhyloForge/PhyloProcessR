#' @title convertAmbiguousConsensus
#'
#' @description Converts IUPAC ambiguity codes in a DNAStringSet alignment to specific
#' unambiguous bases. Each ambiguity code has a corresponding parameter specifying which
#' base to substitute. This is useful for programs that do not accept ambiguous characters.
#' Note: parameter values supplied by the caller are currently overridden inside the
#' function body by fixed defaults (A or T priority).
#'
#' @param alignment a DNAStringSet (or coercible list) containing the sequences to convert.
#'
#' @param convert.r replacement base for the ambiguity code R (A or G). Default c("A", "G").
#'
#' @param convert.y replacement base for the ambiguity code Y (C or T). Default c("C", "T").
#'
#' @param convert.s replacement base for the ambiguity code S (G or C). Default c("G", "C").
#'
#' @param convert.w replacement base for the ambiguity code W (A or T). Default c("A", "T").
#'
#' @param convert.k replacement base for the ambiguity code K (G or T). Default c("G", "T").
#'
#' @param convert.m replacement base for the ambiguity code M (A or C). Default c("A", "C").
#'
#' @param convert.b replacement base for the ambiguity code B (C, G, or T).
#' Default c("C", "G", "T").
#'
#' @param convert.d replacement base for the ambiguity code D (A, G, or T).
#' Default c("A", "G", "T").
#'
#' @param convert.h replacement base for the ambiguity code H (A, C, or T).
#' Default c("A", "C", "T").
#'
#' @param convert.v replacement base for the ambiguity code V (A, C, or G).
#' Default c("A", "C", "G").
#'
#' @return a DNAStringSet with all specified ambiguity codes replaced by their chosen
#' unambiguous bases.
#'
#' @export

convertAmbiguousConsensus = function(alignment = NULL,
                                     convert.r = c("A", "G"),
                                     convert.y = c("C", "T"),
                                     convert.s = c("G", "C"),
                                     convert.w = c("A", "T"),
                                     convert.k = c("G", "T"),
                                     convert.m = c("A", "C"),
                                     convert.b = c("C", "G", "T"),
                                     convert.d = c("A", "G", "T"),
                                     convert.h = c("A", "C", "T"),
                                     convert.v = c("A", "C", "G")) {
  #To do
  # Add randomization
  # Add select best
  #Add column?

  #Debug
 # alignment = non.align
  convert.r = "A"
  convert.y = "T"
  convert.s = "G"
  convert.w = "A"
  convert.k = "T"
  convert.m = "A"
  convert.b = "T"
  convert.d = "T"
  convert.h = "T"
  convert.v = "A"

  #applies across DNASTringSet. Convert formats?
  edit.align = lapply(alignment, function (x) gsub("R|r", convert.r, x))
  edit.align = lapply(edit.align, function (x) gsub("Y|y", convert.y, x))
  edit.align = lapply(edit.align, function (x) gsub("S|s", convert.s, x))
  edit.align = lapply(edit.align, function (x) gsub("W|w", convert.w, x))
  edit.align = lapply(edit.align, function (x) gsub("K|k", convert.k, x))
  edit.align = lapply(edit.align, function (x) gsub("M|m", convert.m, x))
  edit.align = lapply(edit.align, function (x) gsub("B|b", convert.b, x))
  edit.align = lapply(edit.align, function (x) gsub("D|d", convert.d, x))
  edit.align = lapply(edit.align, function (x) gsub("H|h", convert.h, x))
  edit.align = lapply(edit.align, function (x) gsub("V|v", convert.v, x))

  n.align = Biostrings::DNAStringSet(unlist(edit.align))

  return(n.align)

}#end function


