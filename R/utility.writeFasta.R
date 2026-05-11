#' @title writeFasta
#'
#' @description Writes one or more DNA sequences to a FASTA-formatted file.
#'   Sequences can be supplied as a list of character strings or as a single
#'   character string. Long sequences are wrapped at nbchar characters per line.
#'   The file connection can be kept open for appending additional sequences.
#'
#' @param sequences a list of character strings, each element being one
#'   sequence, or a single character string when writing only one sequence.
#'
#' @param names character vector of sequence names corresponding to sequences;
#'   used as FASTA headers (written as ">name").
#'
#' @param file.out path to the output FASTA file to create or append to.
#'
#' @param open character; file connection mode passed to file(). Use "w" to
#'   write (overwrite), "a" to append to an existing file.
#'
#' @param nbchar integer; number of sequence characters per line when wrapping
#'   long sequences. Use a very large value (e.g. 1000000) to write each
#'   sequence on a single line.
#'
#' @param as.string logical; if TRUE each sequence element is treated as a
#'   single character string and converted to a character vector via seqinr::s2c
#'   before writing; if FALSE sequences are expected to already be character
#'   vectors.
#'
#' @return invisibly; writes sequences to file.out in FASTA format.
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

