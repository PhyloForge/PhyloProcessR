#' @title gatherUnlinked
#'
#' @description Assembles a mixed set of gene-level and exon-level alignments for use in
#' unlinked analysis (e.g. coalescent methods). For each gene that has a concatenated gene
#' alignment, that file is copied to the output directory. For genes represented by only a
#' single exon (and therefore absent from the gene alignment directory), the corresponding
#' exon alignment is copied instead. The result is one alignment per locus, avoiding
#' redundancy between gene and exon alignments.
#'
#' @param gene.alignment.directory path to the directory containing gene-level concatenated
#' alignment files (produced by, e.g., \code{concatenateGenes}).
#'
#' @param exon.alignment.directory path to the directory containing individual exon
#' alignment files.
#'
#' @param output.directory path to the directory where the combined set of alignment files
#' will be saved.
#'
#' @param feature.gene.names path to a tab-delimited metadata file with at minimum columns
#' named \code{marker} and \code{gene}, mapping each exon alignment name to a gene name.
#'
#' @param overwrite logical. If TRUE, the output directory is removed and recreated before
#' copying; if FALSE, files are added to an existing directory. Default FALSE.
#'
#' @return Copies alignment files to \code{output.directory}. No value is returned to R.
#'
#' @export

gatherUnlinked = function(gene.alignment.directory = NULL,
                          exon.alignment.directory = NULL,
                          output.directory = NULL,
                          feature.gene.names = NULL,
                          overwrite = FALSE
                          ) {

  #Debug
  # work.dir = "/Volumes/LaCie/data-analysis"
  # setwd(work.dir)
  # gene.alignment.directory = "alignments/untrimmed_genes"
  # exon.alignment.directory = "alignments/untrimmed_all-markers"
  # output.directory = "alignments/untrimmed_all-unique"
  # overwrite = FALSE
  # feature.gene.names = "gene_metadata.txt"

  # Parameter checks
  if (is.null(gene.alignment.directory) == TRUE) {
    stop("Error: a folder of alignments is needed.")
  }
  if (is.null(exon.alignment.directory) == TRUE) {
    stop("Error: a folder of alignments is needed.")
  }

  if (is.null(output.directory) == TRUE) {
    stop("Error: an output file name is needed.")
  }
  if(is.null(feature.gene.names) == TRUE){ stop("Error: a table associating each exon with a gene is needed.") }

  # Check if files exist or not
  if (dir.exists(gene.alignment.directory) == FALSE) {
    return(paste0("Directory of alignments could not be found. Exiting."))
  } # end file check
  
  # Check if files exist or not
  if (dir.exists(exon.alignment.directory) == FALSE) {
    return(paste0("Directory of alignments could not be found. Exiting."))
  } # end file check

  #Checks output overwrite
  if (overwrite == TRUE){
    if (file.exists(paste0(output.directory)) == TRUE) {
      system(paste0("rm -r ", output.directory))
    }
    dir.create(output.directory)
  } else {
    dir.create(output.directory)
  }#end overwrite if

  # Gets list of alignments
  gene.files = list.files(gene.alignment.directory, full.names = FALSE, recursive = TRUE)
  exon.files = list.files(exon.alignment.directory, full.names = FALSE, recursive = TRUE)
  exon.data = data.table::fread(file = feature.gene.names, header = TRUE)
  
  single.data = exon.data[!exon.data$gene %in% gsub("\\..*", "", gene.files),]
  
  #Copies the genes over
  for (i in seq_along(gene.files)){
    system(paste0("cp ", gene.alignment.directory, "/", gene.files[i], " ", output.directory, "/", gene.files[i]))
  }
  
  #Copies the remaining exons over
  exon.copy = exon.files[gsub("\\..*", "", exon.files) %in% single.data$marker]

  # Copies the genes over
  for (i in seq_along(exon.copy)) {
    system(paste0("cp ", exon.alignment.directory, "/", exon.copy[i], " ", output.directory, "/", exon.copy[i]))
  }

}#end function