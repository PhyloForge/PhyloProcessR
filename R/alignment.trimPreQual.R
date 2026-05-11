#' @title trimPreQual
#'
#' @description Wrapper function for running HmmCleaner (via HmmCleaner.pl) on a single alignment. The alignment is written to a temporary fasta file, HmmCleaner is called with the selected options, and the cleaned alignment is read back and passed through trimAlignmentColumns to remove any fully-gap columns introduced by HmmCleaner. If HmmCleaner fails or produces no output the original alignment is returned unchanged. HmmCleaner must be installed and accessible.
#'
#' @param alignment a DNAStringSet containing the aligned sequences to clean
#'
#' @param specificity if TRUE, run HmmCleaner with the --specificity flag for more sensitive but slower cleaning
#'
#' @param large if TRUE and the alignment has 50 or more sequences, run HmmCleaner with the --large flag for faster processing
#'
#' @param hmmcleaner.path path to the HmmCleaner.pl script or the directory containing it; defaults to "prequal"
#'
#' @param quiet if TRUE, suppress HmmCleaner screen output
#'
#' @param delete.temp if TRUE, delete the temporary fasta files created during the HmmCleaner run
#'
#' @return a DNAStringSet of the HmmCleaner-cleaned alignment with all-gap columns removed; returns the original alignment unchanged if HmmCleaner produces no output
#'
#' @export

trimPreQual = function(alignment = NULL,
                       specificity = TRUE,
                       large = FALSE,
                       hmmcleaner.path = "prequal",
                       quiet = FALSE,
                       delete.temp = TRUE){

  #Debug section
  alignment = non.align
  specificity = TRUE
  large = FALSE
  hmmcleaner.path = hmm.path
  quiet = FALSE
  delete.temp = TRUE

  #Creates random name and saves it
  write.align = as.list(as.character(alignment))
  input.file = paste0("temp", sample(1:10000, 1), ".fa")
  writeFasta(sequences = write.align,
             names = names(write.align),
             file.out = input.file,
             nbchar = 1000000,
             as.string = T)

  #For datasets with 50 or more samples
  if (large == TRUE & length(alignment) >= 50 & specificity == FALSE){
    system(paste0("HmmCleaner.pl ", input.file, " --large"),
           ignore.stdout = quiet, ignore.stderr = quiet)
  }# end if

  #Specificity function
  if (specificity == TRUE & large == FALSE){
    system(paste0("HmmCleaner.pl ", input.file, " --specificity"),
           ignore.stdout = quiet, ignore.stderr = quiet)
  }#end if

  #Specificity function
  if (large == TRUE & length(alignment) >= 50 & specificity == TRUE){
    system(paste0("HmmCleaner.pl ", input.file, " --large --specificity"),
           ignore.stdout = quiet, ignore.stderr = quiet)
  }#end if

  if (file.exists(paste0(gsub(".fa$", "", input.file), "_hmm.fasta")) == TRUE){
    hmm.align = Rsamtools::scanFa(Rsamtools::FaFile(paste0(gsub(".fa$", "", input.file), "_hmm.fasta")))
  } else {
    if (delete.temp == TRUE){ system(paste0("rm ", gsub(".fa$", "", input.file), "*")) }
    return(alignment)
  }#end file check

  if (delete.temp == TRUE){ system(paste0("rm ", gsub(".fa$", "", input.file), "*")) }

  #Removes the edge gaps
  out.align = trimAlignmentColumns(alignment = hmm.align,
                                   min.gap.percent = 100)

  # #Remove gap only alignments
  # gap.align = strsplit(as.character(out.align), "")
  # gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
  # gap.rem = gap.count[gap.count == Biostrings::width(out.align)[1]]

  return(out.align)

}#end function
