#' @title trimHmmCleaner
#'
#' @description wrapper function for running HmmCleaner. Must be installed.
#'
#' @param alignment alignment in DNAStringSet format
#'
#' @param specificity TRUE to enable to more sensitive but slower "specificity" hmmCleaner option
#'
#' @param large TRUE to enable the faster "large" dataset option for large datasets. Will not run "large" if alignment has less than 50 sequences
#'
#' @param hmmcleaner.path Absolute path to hmmCleaner if R cannot find it in your path
#'
#' @param quiet TRUE to supress HmmCleaner screen output
#'
#' @param delete.temp TRUE to delete temporary files made by HmmCleaner
#'
#' @return returns a data.table with the raw summary statistics calculated for each alignment in the alignment.path. A csv file can optionally be saved by giving a file name to file.export
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

trimHmmCleaner = function(alignment = NULL,
                          specificity = TRUE,
                          large = FALSE,
                          hmmcleaner.path = "HmmCleaner.pl",
                          quiet = FALSE,
                          delete.temp = TRUE){

  #Debug section
  # alignment = non.align
  # specificity = TRUE
  # large = FALSE
  # hmmcleaner.path = hmm.path
  # quiet = FALSE
  # delete.temp = TRUE

  if (length(alignment) <= 3){ return(alignment) }

  #Creates random name and saves it
  write.align = as.list(as.character(alignment))
  input.file = paste0("temp", sample(1:1000000, 1), ".fa")
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
           ignore.stdout = quiet, ignore.stderr = F)
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
