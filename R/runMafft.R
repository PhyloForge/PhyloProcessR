#' @title runMafft
#'
#' @description Function for gather summary statistics on your alignment. Can be used for filtering or summarizing data.
#'
#' @param sequence.data path to a folder of sequence alignments in phylip format.
#'
#' @param add.contigs contigs are added into existing alignment if algorithm is "add"
#'
#' @param algorithm algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param adjust.direction TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param save.name if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param cleanup.files give a save name if you wnat to save the summary to file.
#'
#' @param quiet TRUE to supress mafft screen output

#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
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

runMafft = function(sequence.data = NULL,
                    add.contigs = NULL,
                    algorithm = "localpair",
                    adjust.direction = TRUE,
                    save.name = NULL,
                    threads = 1,
                    cleanup.files = TRUE,
                    quiet = TRUE){

  #save.name<-locus.save.name
  #algorithm = "localpair"
  # unaligned.contigs<-intron.align
  save.contigs = as.list(as.character(sequence.data))
  if (is.null(save.name) == T) { save.name = paste(sample(LETTERS, 5, replace = T), collapse = "")}
  if (adjust.direction == T){ adjust.direction = "--adjustdirection" } else { adjust.direction = "" }

  #Adds a sequence into the alignment. Saves much computation.
  if (algorithm == "add"){
    #Saves to folder to run with mafft
    writeFasta(sequences = save.contigs, names = names(save.contigs),
               paste0(save.name, ".fa"), nbchar = 1000000, as.string = T)

    #Saves to folder to run with mafft
    add.save = as.list(as.character(add.contigs))
    writeFasta(sequences = add.save, names = names(add.save),
               "add_sequences.fa", nbchar = 1000000, as.string = T)

    #Runs MAFFT to align
    system(paste0("mafft --",algorithm, " add_sequences.fa --maxiterate 1000 ", save.name, ".fa > ",
                  save.name, "_align.fa"), ignore.stderr = T)

    alignment = Rsamtools::scanFa(Rsamtools::FaFile(paste(save.name, "_align.fa")))   # loads up fasta file
    unlink(paste0(save.name, ".fa"))
    unlink("add_sequences.fa")

  }#end -add

  #Does Regular MAFFT Local Pair
  if (algorithm == "localpair"){
    #Saves to folder to run with mafft
    writeFasta(sequences = save.contigs, names = names(save.contigs),
               paste(save.name, ".fa", sep = ""), nbchar = 1000000, as.string = T)

    #Runs MAFFT to align
    system(paste0("mafft --",algorithm, " --maxiterate 1000 ", adjust.direction, " --quiet --op 3 --ep 0.123",
                  " --thread ", threads, " ", save.name, ".fa > ", save.name, "_align.fa"), ignore.stderr = T)

    alignment = Rsamtools::scanFa(Rsamtools::FaFile(paste0(save.name, "_align.fa")))   # loads up fasta file
    unlink(paste0(save.name, ".fa"))
  }#end local pair

  if (cleanup.files == T){
    unlink(paste0(save.name, "_align.fa"))
    return(alignment)
  } else { return(alignment) }
}#function end


