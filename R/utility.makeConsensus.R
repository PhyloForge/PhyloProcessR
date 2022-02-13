#' @title makeConsensus
#'
#' @description Function for making a consensus sequence from alignment
#'
#' @param alignment alignment in DNAStringSet format
#'
#' @param method method to apply; majority = majority consensus, threshold = consensus from threshold percent of column; IUPAC and profile
#'
#' @param threshold threshold to use for "threshold" method
#'
#' @param warn.non.IUPAC warn if characters are not IUPAC
#'
#' @param type DNA or RNA sequence data
#'
#' @return returns a data.table with the raw summary statistics calculated for each alignment in the alignment.path. A csv file can optionally be saved by giving a file name to file.export
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

makeConsensus = function(alignment = NULL,
                         method = c("majority", "threshold", "IUPAC", "profile"),
                         threshold = 0.6,
                         warn.non.IUPAC = FALSE,
                         remove.gaps = TRUE,
                         type = c("DNA", "RNA")) {

  #input.alignment<-temp.align
  #Converts alignment to matrix of characters to be used
  new.align<-strsplit(as.character(alignment), "")
  align.in<-matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)

  #Does based on method
  method <- match.arg(method)

  if (method == "IUPAC") {
    type <- match.arg(type)
    res <- apply(align.in, 2, bma, warn.non.IUPAC = warn.non.IUPAC,
                 type = type)
    names(res) <- NULL
  }
  if (method == "majority") {
    majority <- function(x) names(which.max(table(x)))
    res <- apply(align.in, 2, majority)
    names(res) <- NULL
  }
  if (method == "profile") {
    obsvalue <- levels(factor(align.in))
    nrow <- length(obsvalue)
    row.names(align.in) <- NULL
    res <- apply(align.in, 2, function(x) table(factor(x, levels = obsvalue)))
  }
  if (method == "threshold") {
    profile <- consensus(align.in, method = "profile")
    profile.rf <- apply(profile, 2, function(x) x/sum(x))
    res <- rownames(profile.rf)[apply(profile.rf, 2, which.max)]
    res <- ifelse(apply(profile.rf, 2, max) >= threshold,
                  res, NA)
    names(res) <- NULL
  }

  out.consensus<-Biostrings::DNAStringSet(paste0(res, collapse = ""))
  if (remove.gaps == TRUE){
    out.consensus = Biostrings::DNAStringSet(gsub("-", "", out.consensus))
  }#end if

  names(out.consensus)<-"Consensus_Sequence"
  return(out.consensus)
}
