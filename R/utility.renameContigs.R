#' @title renameContigs
#'
#' @description Renames contig FASTA files by copying them from a source
#'   directory to a new directory with names taken from a rename spreadsheet.
#'   The spreadsheet maps current file names to desired sample names. The
#'   function supports custom column header names for the source and destination
#'   name columns.
#'
#' @param rename.spreadsheet path to a CSV file mapping existing contig file
#'   names to new sample names. By default the first column is used as the
#'   source name and the second as the new name; this is controlled by headers.
#'
#' @param contig.directory path to the directory containing the input contig
#'   FASTA files (.fa or .fasta).
#'
#' @param output.directory path to the directory where renamed contig files
#'   will be saved.
#'
#' @param skip.not.found logical; if TRUE samples in the spreadsheet that cannot
#'   be located in contig.directory are silently skipped. If FALSE an error is
#'   raised.
#'
#' @param headers character vector of length 2 specifying the column names in
#'   rename.spreadsheet to use as source and destination names respectively
#'   (e.g. c("File", "Sample")). NULL uses the first and second columns.
#'
#' @param overwrite logical; if TRUE the output directory must not exist (the
#'   function raises an error if it does and overwrite is FALSE); if TRUE the
#'   directory is deleted and recreated.
#'
#' @return invisibly; copies renamed contig FASTA files to output.directory.
#'
#' @export

renameContigs = function(rename.spreadsheet = NULL,
                         contig.directory = NULL,
                         output.directory = NULL,
                         skip.not.found = FALSE,
                         headers = NULL,
                         overwrite = FALSE){

  # ### Example usage
  # rename.spreadsheet = "/Users/chutter/Dropbox/SharewithCarl/file_rename.csv"
  # output.directory = "/Users/chutter/Dropbox/SharewithCarl/renamed_contigs"
  # contig.directory = "/Users/chutter/Dropbox/SharewithCarl/draft-assemblies"
  # overwrite = FALSE
  # headers = c("File", "Sample")
  # skip.not.found = FALSE

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }

    if (overwrite == FALSE){ stop("Overwrite == FALSE and directory exists.") }
  }#end else

  sample.data = read.csv(rename.spreadsheet)
  all.contigs = list.files(contig.directory)
  all.contigs = all.contigs[grep("fasta$|fa$", all.contigs)]

  if (nrow(sample.data) == 0){ return("no contigs to rename in spreadsheet.") }
  if (length(all.contigs) == 0){ return("no contigs to rename in directory.") }

  for (i in 1:nrow(sample.data)){

    if (is.null(headers) == TRUE){
      sample.contigs = all.contigs[grep(as.character(sample.data[i,1]), gsub(".fa$|.fasta$", "", all.contigs))]
    }

    if (is.null(headers) == FALSE){
      sample.contigs = all.contigs[grep(sample.data[,colnames(sample.data) == headers[1]][i], gsub(".fa$|.fasta$", "", all.contigs))]
    }

    if (length(sample.contigs) == 0){
      if (skip.not.found == FALSE) { stop("Sample in spreadsheet not found in contigs.") }
      if (skip.not.found == TRUE) { next }
    }

    #Save the read files with the new names in the new directory
    if (is.null(headers) == TRUE){
      new.name = sample.data[grep(gsub(".fa$|.fasta$", "", sample.contigs), as.character(sample.data[,1])),][2]
    }

    if (is.null(headers) == FALSE){
      new.name = sample.data[grep(gsub(".fa$|.fasta$", "", sample.contigs), sample.data[,colnames(sample.data) == headers[1]]),]
      new.name = new.name[,colnames(sample.data) == headers[2]]
    }


    system(paste0("cp ", contig.directory, "/", sample.contigs, " ",
                  paste0(output.directory, "/",  new.name, ".fa")))

  }#end i loop

}#end function






# END SCRIPT


