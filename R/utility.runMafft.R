#' @title runMafft
#'
#' @description Aligns sequences using MAFFT. Supports local-pair, global-pair,
#'   and "add" modes. In "add" mode, new sequences (add.contigs) are inserted
#'   into an existing alignment (sequence.data) without re-aligning all
#'   sequences. All other modes perform a full multiple sequence alignment of
#'   sequence.data. Temporary FASTA files are written and optionally cleaned up
#'   after alignment.
#'
#' @param sequence.data a DNAStringSet of sequences to align (or the existing
#'   alignment when algorithm = "add").
#'
#' @param add.contigs a DNAStringSet of sequences to add to the alignment when
#'   algorithm = "add"; ignored for other algorithms.
#'
#' @param algorithm character; the MAFFT alignment strategy to use. "localpair"
#'   and "globalpair" perform a full alignment of sequence.data. "add" adds
#'   add.contigs to sequence.data without realigning existing sequences.
#'
#' @param dna.type character; sequence type passed to MAFFT (e.g. "nuc" for
#'   nucleotide, "amino" for protein).
#'
#' @param adjust.direction logical; if TRUE passes --adjustdirection to MAFFT
#'   to automatically reverse-complement sequences that align better in the
#'   reverse orientation.
#'
#' @param save.name base name used for the temporary FASTA files written for
#'   MAFFT; if NULL a random five-letter string is generated.
#'
#' @param threads integer number of CPU threads to pass to MAFFT (used for
#'   non-"add" algorithms).
#'
#' @param cleanup.files logical; if TRUE the temporary FASTA file produced by
#'   MAFFT is deleted after reading the alignment.
#'
#' @param mafft.path system path to the directory containing the mafft
#'   executable; NULL searches the system PATH.
#'
#' @param quiet logical; if TRUE MAFFT stderr is suppressed.
#'
#' @return a DNAStringSet containing the aligned sequences.
#'
#' @export

runMafft = function(sequence.data = NULL,
                    add.contigs = NULL,
                    algorithm = c("localpair", "globalpair", "add"),
                    dna.type = "nuc",
                    adjust.direction = TRUE,
                    save.name = NULL,
                    threads = 1,
                    cleanup.files = TRUE,
                    mafft.path = NULL,
                    quiet = TRUE){

  #save.name<-locus.save.name
  #algorithm = "localpair"
  # unaligned.contigs<-intron.align

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  save.contigs = as.list(as.character(sequence.data))
  if (is.null(save.name) == T) { save.name = paste(sample(LETTERS, 5, replace = T), collapse = "")}
  if (adjust.direction == T){ adjust.direction = "--adjustdirection " } else { adjust.direction = "" }

  #Adds a sequence into the alignment. Saves much computation.
  if (algorithm == "add"){
    #Saves to folder to run with mafft
    writeFasta(sequences = save.contigs, names = names(save.contigs),
               paste0(save.name, ".fa"), nbchar = 1000000, as.string = T)

    #Saves to folder to run with mafft
    add.save = as.list(as.character(add.contigs))
    writeFasta(sequences = add.save, names = names(add.save),
               paste0(save.name, "_add_sequences.fa"), nbchar = 1000000, as.string = T)

    #Runs MAFFT to align
    system(paste0(mafft.path, "mafft --",algorithm, " ", save.name,
                  "_add_sequences.fa --maxiterate 1000 --", dna.type, " ",
                  adjust.direction, save.name, ".fa > ",
                  save.name, "_align.fa"), ignore.stderr = quiet)

    alignment = Biostrings::readDNAStringSet(paste0(save.name, "_align.fa"))   # loads up fasta file
    unlink(paste0(save.name, ".fa"))
    unlink(paste0(save.name, "_add_sequences.fa"))

  }#end -add

  #Does Regular MAFFT Local Pair
  if (algorithm != "add"){
    #Saves to folder to run with mafft
    writeFasta(sequences = save.contigs, names = names(save.contigs),
               paste(save.name, ".fa", sep = ""), nbchar = 1000000, as.string = T)

    #Runs MAFFT to align
    system(paste0(mafft.path, "mafft --",algorithm, " --maxiterate 1000 --", dna.type, " ",
                  adjust.direction, "--quiet --op 3 --ep 0.123",
                  " --thread ", threads, " ", save.name, ".fa > ", save.name, "_align.fa"),
           ignore.stderr = quiet)

    alignment = Biostrings::readDNAStringSet(paste0(save.name, "_align.fa"))   # loads up fasta file
    unlink(paste0(save.name, ".fa"))
  }#end local pair

  if (cleanup.files == T){
    unlink(paste0(save.name, "_align.fa"))
    return(alignment)
  } else { return(alignment) }
}#function end


