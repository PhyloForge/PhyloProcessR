#' @title makeConsensus
#'
#' @description Derives a consensus sequence from a multiple sequence alignment
#'   stored as a DNAStringSet. Four methods are available: majority rule,
#'   threshold-based, IUPAC ambiguity codes, or a positional frequency profile.
#'   Gap characters are optionally removed from the consensus.
#'
#' @param alignment a DNAStringSet containing an aligned set of sequences (all
#'   must have equal length).
#'
#' @param method character; the consensus algorithm to apply. "majority" assigns
#'   the most frequent base at each column; "threshold" assigns a base only when
#'   it exceeds the specified frequency threshold; "IUPAC" uses seqinr::bma to
#'   produce IUPAC ambiguity codes; "profile" returns a positional frequency
#'   matrix rather than a sequence.
#'
#' @param threshold numeric (0-1); the minimum column frequency required to call
#'   a base when method = "threshold".
#'
#' @param warn.non.IUPAC logical; if TRUE a warning is issued for non-IUPAC
#'   characters (only used when method = "IUPAC").
#'
#' @param remove.gaps logical; if TRUE gap characters ("-") are removed from
#'   the final consensus sequence.
#'
#' @param type character; "DNA" or "RNA", passed to seqinr::bma when method =
#'   "IUPAC".
#'
#' @return a DNAStringSet of length 1 named "Consensus_Sequence" containing the
#'   consensus sequence. When method = "profile" a frequency matrix is returned
#'   instead.
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
    res <- apply(align.in, 2, seqinr::bma, warn.non.IUPAC = warn.non.IUPAC,
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
