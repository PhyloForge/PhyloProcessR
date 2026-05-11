#' @title shortenFastaHeaders
#'
#' @description Truncates FASTA sequence headers to a maximum number of
#'   characters for all FASTA files in a directory. This is useful when
#'   downstream tools (e.g. some alignment or tree programs) impose a header
#'   length limit. Each input file is read, headers are trimmed with
#'   strtrim(), and the result is saved to the output directory under the same
#'   file name.
#'
#' @param fasta.directory path to a directory containing FASTA files whose
#'   headers should be shortened.
#'
#' @param output.directory path to the directory where files with shortened
#'   headers will be saved.
#'
#' @param number.characters integer maximum number of characters to retain in
#'   each sequence header (headers longer than this are truncated).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before writing.
#'
#' @param quiet logical; currently unused, reserved for future output control.
#'
#' @return invisibly; writes FASTA files with shortened headers to
#'   output.directory.
#'
#' @export

shortenFastaHeaders = function(fasta.directory = NULL,
                               output.directory = "shortened-headers",
                               number.characters = 50,
                               overwrite = FALSE,
                               quiet = TRUE) {


  #Read in basic genome info
  #fasta.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/proteomes"
  #output.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Mantellidae_Subfamily/Mantella_Genome/maker-analysis/shortened-headers"
  #number.characters = 70

  if (is.null(fasta.directory) == T){ stop("A directory of genome(s) is needed.") }

  #Overwrites
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

