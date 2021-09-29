#' @title writeFasta
#'
#' @description Function for writing sequences or alignment in R to fasta format
#'
#' @param sequences the sequences in some format to do
#'
#' @param names names from the sequences
#'
#' @param file.out name of the file to be saved
#'
#' @param open keep the connection open = a to append, or 'w" to write once
#'
#' @param nbchar number of characters per line if interleaving
#'
#' @param as.string save as a string or not
#'
#' @return saves the alignment as a fasta file
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

##### saves alignment as a phylip file
writeFasta = function(sequences,
                      names,
                      file.out,
                      open = "w",
                      nbchar = 60,
                      as.string = FALSE){
  #
  # Open output file:
  #
  outfile <- file(description = file.out, open = open)

  #
  # Function to write one sequence in output file:
  #
  write.oneseq <- function(sequence, name, nbchar, as.string){
    writeLines(paste(">", name, sep = ""), outfile)
    if(as.string) sequence <- seqinr::s2c(sequence)
    l <- length(sequence)
    q <- floor(l/nbchar)
    r <- l - nbchar*q
    if(q > 0){
      sapply(seq_len(q), function(x) writeLines(seqinr::c2s(sequence[(nbchar*(x - 1) + 1):(nbchar*x)]), outfile))
    }
    if(r > 0){
      writeLines(seqinr::c2s(sequence[(nbchar*q + 1):l]), outfile)
    }
  }

  #
  # Write all sequences in output file:
  #
  if(!is.list(sequences)){
    write.oneseq(sequence = sequences, name = names, nbchar = nbchar, as.string = as.string)
  } else {
    n.seq <- length(sequences)
    sapply(seq_len(n.seq), function(x) write.oneseq(sequence = as.character(sequences[[x]]),
                                                    name = names[x], nbchar = nbchar, as.string = as.string))
  }
  #
  # Close output file:
  #
  close(outfile)
}

