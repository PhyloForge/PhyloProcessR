#' @title gatherUnlinked
#'
#' @description Function for concatenating a large number of alignments
#'
#' @param alignment.folder folder that contains aligmnents to be concatenated
#'
#' @param output.folder output file name
#'
#' @param exon.gene.names output file name
#'
#' @param input.format input file format. Save three types: phylip, nexus, and fasta
#'
#' @param output.format output file format. Save three types: phylip, nexus, and fasta
#'
#' @param remove.reverse TRUE to remove "_R_" placed before reversed sequences in some alignments. Default FALSE.
#'
#' @param overwrite TRUE to overwrite file. Default FALSE.
#'
#' @param threads TRUE to overwrite file. Default FALSE.
#'
#' @param memory TRUE to overwrite file. Default FALSE.
#'
#' @return saves to file concatenated alignments and partition files delimiting the coordinates of each indidividual marker
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
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
      system(paste0("rm -r ", output.folder))
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