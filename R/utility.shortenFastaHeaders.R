#' @title shortenFastaHeaders
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

shortenFastaHeaders = function(fasta.directory = NULL,
                               output.directory = "shortened-headers",
                               number.characters = 70,
                               overwrite = FALSE,
                               resume = TRUE,
                               quiet = TRUE) {


  #Read in basic genome info
  #fasta.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/proteomes"
  #output.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/shortened-headers"
  #number.characters = 70

  if (is.null(fasta.directory) == T){ stop("A directory of genome(s) is needed.") }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  fasta.files = list.files(fasta.directory)

  for (i in 1:length(fasta.files)){

    old.file = Biostrings::readAAStringSet(paste0(fasta.directory, "/", fasta.files[i]), format = "fasta")
    names(old.file) = strtrim(x = names(old.file), width = number.characters)

    #Creates random name and saves it
    write.file = as.list(as.character(old.file))
    writeFasta(sequences = write.file,
               names = names(write.file),
               file.out = paste0(output.directory, "/", fasta.files[i]),
               nbchar = 100000000,
               as.string = T)

  }# end i loop

}#end funtion

