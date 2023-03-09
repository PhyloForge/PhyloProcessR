#' @title renameContigs
#'
#' @description Function for downloading data off your dropbox account
#'
#' @param sample.spreadsheet spreadsheet path or table of file names to download
#'
#' @param dropbox.directory the directory that contains these target files
#'
#' @param dropbox.token your dropbox token path, created the first time you run the function but needs interactivity for log-in. Can be done once on your computer and the token can be moved around.
#'
#' @param output.directory the directory to save the files you are downloading
#'
#' @param skip.not.found the directory to save the files you are downloading
#'
#' @param overwrite whether to overwrite or not
#'
#' @return saves dropbox files to output directory
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
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


