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
    stop("Error: gene.alignment.directory is required.")
  }
  if (is.null(exon.alignment.directory) == TRUE) {
    stop("Error: exon.alignment.directory is required.")
  }

  if (is.null(output.directory) == TRUE) {
    stop("Error: an output file name is needed.")
  }
  if(is.null(feature.gene.names) == TRUE){ stop("Error: a table associating each exon with a gene is needed.") }

  # Check if files exist or not
  if (dir.exists(gene.alignment.directory) == FALSE) {
    stop("gene.alignment.directory not found. Please check the path.")
  }
  if (dir.exists(exon.alignment.directory) == FALSE) {
    stop("exon.alignment.directory not found. Please check the path.")
  }

  exon.data = .readGeneMetadata(feature.gene.names)
  gene.files = .alignmentFiles(gene.alignment.directory)
  exon.files = .alignmentFiles(exon.alignment.directory)
  if (length(gene.files) == 0 && length(exon.files) == 0) {
    stop("No alignment files were found in the gene or exon directory.")
  }
  matched.markers = .alignmentId(exon.files) %in% exon.data$marker
  if (length(exon.files) > 0 && !any(matched.markers) && length(gene.files) == 0) {
    stop("No alignment names match the gene metadata marker column.")
  }

  # Checks output overwrite after validating all inputs.
  if (overwrite == TRUE){
    if (dir.exists(output.directory) == TRUE) unlink(output.directory, recursive = TRUE)
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (!dir.exists(output.directory)) dir.create(output.directory, recursive = TRUE)
  }#end overwrite if

  single.data = exon.data[!gene %in% .alignmentId(gene.files)]
  
  #Copies the genes over
  for (i in seq_along(gene.files)){
    dest = paste0(output.directory, "/", gene.files[i])
    if (overwrite == FALSE && file.exists(dest)) { next }
    .copyAlignment(file.path(gene.alignment.directory, gene.files[i]), dest,
                   overwrite = overwrite)
  }

  #Copies the remaining single-exon loci over
  exon.copy = exon.files[.alignmentId(exon.files) %in% single.data$marker]

  for (i in seq_along(exon.copy)) {
    dest = paste0(output.directory, "/", exon.copy[i])
    if (overwrite == FALSE && file.exists(dest)) { next }
    .copyAlignment(file.path(exon.alignment.directory, exon.copy[i]), dest,
                   overwrite = overwrite)
  }

}#end function
