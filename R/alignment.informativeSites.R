#' @title informativeSites
#'
#' @description Counts or calculates the proportion of parsimony-informative sites in an
#' alignment. A site is parsimony-informative if at least two different character states
#' each appear in at least two sequences, excluding gap (\code{-}), missing (\code{?}),
#' and N characters (and optionally IUPAC ambiguity codes). The function operates on an
#' alignment in DNAbin matrix format.
#'
#' @param alignment alignment in ape DNAbin matrix format (rows = taxa, columns = sites).
#'
#' @param count logical. If TRUE (default), returns the integer count of parsimony-
#' informative sites. If FALSE, returns the proportion of informative sites relative to
#' alignment length.
#'
#' @param variable.sites logical. Currently passed to the function but the implementation
#' always computes parsimony-informative sites. Included for future extension. Default FALSE.
#'
#' @param ambiguities logical. If TRUE (default), IUPAC ambiguity codes are treated as
#' valid characters. If FALSE, ambiguity codes are added to the exclusion list and not
#' counted.
#'
#' @return an integer count (if \code{count = TRUE}) or a numeric proportion rounded to
#' three decimal places (if \code{count = FALSE}) of parsimony-informative sites.
#'
#' @examples
#' \dontrun{
#' n.pis <- informativeSites(alignment = as.matrix(ape::as.DNAbin(my.align)))
#' }
#'
#' @export

#### MAKE NEW FUNCTION THAT OUTPUTS ALL. UPDATE OTHER FUNCTIONS THAT USE IT

#Calculates informative sites
informativeSites = function(alignment = NULL,
                            count = TRUE,
                            variable.sites = FALSE,
                            ambiguities = TRUE) {

  #Helper function to use with apply
  column.pars = function(x) {
    x = table(x)
    x = x[x > 1]
    if (length(x[!names(x) %in% n]) > 1){ return(TRUE) } else { return(FALSE) }
  }#end function

  #characters to exclude
  n = c("-", "?", "n")
  if (ambiguities == FALSE){
    n = append(n, c( "r", "y", "k", "m", "s", "w", "b", "d", "h", "v")) }

  #alignment length
  x.len = dim(alignment)[2]
  #goes through each column and sees if they are different
  alignment = as.character(alignment)
  col.pis = apply(alignment, 2, column.pars)
  out = length(col.pis[col.pis == TRUE])

  if (count != TRUE){ out = round(out/x.len, digits = 3) }
  return(out)

}#end informative sites function
