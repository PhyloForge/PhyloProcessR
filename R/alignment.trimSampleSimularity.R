#' @title trimSampleSimilarity
#'
#' @description Removes entire samples from an alignment that are too divergent from the majority-rule consensus sequence. A consensus is computed from the alignment, pairwise distances from the consensus are calculated for each sample, and any sample whose distance meets or exceeds the similarity.threshold is removed. If any samples are removed and realign.mafft is TRUE, the remaining sequences are realigned with MAFFT (localpair). Alignments with two or fewer sequences are returned unmodified.
#'
#' @param alignment a DNAStringSet containing the aligned sequences to filter
#'
#' @param similarity.threshold pairwise distance threshold (0-1): samples with a distance from the consensus at or above this value are removed
#'
#' @param realign.mafft if TRUE and samples were removed, realign the remaining sequences with MAFFT (localpair algorithm) before returning
#'
#' @param mafft.path system path to the MAFFT executable directory; NULL to use the system PATH
#'
#' @return a DNAStringSet with overly divergent samples removed, optionally realigned with MAFFT
#'
#' @export

trimSampleSimilarity = function(alignment = NULL,
                                similarity.threshold = 0.4,
                                realign.mafft = TRUE,
                                mafft.path = NULL) {

  #Debug section
  # alignment = non.align
  # similarity.threshold = 0.45
  # realign.mafft = TRUE
  # mafft.path  = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"

  if (length(alignment) <= 2){ return(alignment) }

  #Make a consensus sequence
  con.seq = makeConsensus(alignment = alignment,
                          method = "majority",
                          type = "DNA",
                          remove.gaps = FALSE)

  #Finds the alignemnt pairwise distance from the target
  diff = pairwiseDistanceTarget(alignment = alignment,
                                target = con.seq)

  #Gets the divergence to make sure not crazy
  bad.seqs = names(diff)[which(diff >= similarity.threshold)]
  rem.align = alignment[!names(alignment) %in% bad.seqs]

  # Moves onto next loop in there are no good sequences
  if (length(rem.align) <= 2){ return(rem.align) }

  ### realign if bad seqs removed
  if (length(bad.seqs) != 0 & realign.mafft == TRUE){
    #Runs mafft
    rem.align = runMafft(sequence.data = rem.align,
                         add.contigs = NULL,
                         algorithm = "localpair",
                         adjust.direction = TRUE,
                         save.name = NULL,
                         threads = 1,
                         cleanup.files = TRUE,
                         quiet = TRUE,
                         mafft.path = mafft.path)

    names(rem.align) = gsub(pattern = "^_R_", replacement = "", x = names(rem.align))

  } # end bad.seqs if

  return(rem.align)

}#end function


