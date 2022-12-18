#' @title makeUniqueFastaHeaders
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param genome.directory path to a folder of sequence alignments in phylip format.
#'
#' @param output.directory available input alignment formats: fasta or phylip
#'
#' @param threads contigs are added into existing alignment if algorithm is "add"
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite TRUE to supress mafft screen output
#'
#' @param resume TRUE to supress mafft screen output
#'
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

makeUniqueFastaHeaders = function(fasta.file = NULL,
                                  output.name = "unique-headers",
                                  type = c("append-number", "rename"),
                                  overwrite = FALSE,
                                  quiet = TRUE) {


  #Read in basic genome info
  #fasta.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/proteomes"
  #output.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/shortened-headers"
  #number.characters = 70
  #fasta.file = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/venom_loci_updated_Mar12_cdhit95.fa"
  #output.name = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/venom_loci_updated_Mar12_cdhit95_unique.fa"
  #type = "append-number"

  if (is.null(fasta.file) == T){ stop("A directory of genome(s) is needed.") }

  #Overwrites
  if (file.exists(output.name) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm ", output.directory))
    } else { stop("file exists and overwrite = FALSE") }
  }

  old.file = Biostrings::readAAStringSet(fasta.file, format = "fasta")

  dup.entry = old.file[duplicated(names(old.file)) == T]
  keep.entry = old.file[duplicated(names(old.file)) == F]

  if (type == "append-number"){
    names(dup.entry) = paste0(names(dup.entry), "_", seq(1:length(dup.entry)))
  }

  save.entry = append(keep.entry, dup.entry)
  save.entry = save.entry[order(names(save.entry))]

  #Renames them all Sequence_X
  if (type == "rename"){
    names(save.entry) = paste0("Sequence_", stringr::str_pad(string = seq(1:length(save.entry)),
                                                             width = nchar(length(save.entry)),
                                                             side = "left",
                                                             pad = "0") )
  }

  #Creates random name and saves it
  write.file = as.list(as.character(save.entry))
  writeFasta(sequences = write.file,
             names = names(write.file),
             file.out = output.name,
             nbchar = 100000000,
             as.string = T)

}#end funtion


