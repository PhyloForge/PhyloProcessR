#' @title trimSampleSimilarity
#'
#' @description Function for trimming out divergent sample sequence segments in an alignment
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param similarity.threshold divergence threshold from consensus to remove a sample
#'
#' @param realign.mafft TRUE to realign alignment with mafft after removing samples
#'
#' @return returns DNAStringSet of trimmed alignment
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
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

    names(rem.align) = gsub(pattern = "_R_", replacement = "", x = names(rem.align))

  } # end bad.seqs if

  return(rem.align)

}#end function


