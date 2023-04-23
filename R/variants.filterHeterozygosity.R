#' @title filterHeterozygosity
#'
#' @description Calculates the number or proportion of parsimony informative sites in an alignment
#'
#' @param alignment alignment in ape DNABin or a matrix format
#'
#' @param count Whethe to return the count of parsimoney informative sites (TRUE) or the proportion (FALSE)
#'
#' @param ambiguities Whether to consider ambiguities (TRUE) or not (FALSE)
#'
#' @return plots the phylogenetic tree and selected data associated with an AstralPlane object. Can optionally be saved to file as a PDF by giving save.file a file name.
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


#Calculates informative sites
filterHeterozygosity = function(iupac.directory = NULL,
                                output.directory = NULL,
                                removed.directory = NULL,
                                threshold = 0.05,
                                overwrite = FALSE
                                ) {

  setwd("/Volumes/LaCie/Mantellidae")
  iupac.directory = "data-analysis/variant-calling/iupac-contigs"
  removed.directory = "data-analysis/contigs/high-het-contigs"
  output.directory = "data-analysis/contigs/filtered-contigs"
  threshold = 0.05
  overwrite = FALSE
  gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  # Initial checks
  if (iupac.directory == output.directory) {
    stop("You should not overwrite the original contigs.")
  }

  if (iupac.directory == removed.directory) {
    stop("You should not overwrite the original contigs.")
  }

  if (is.null(output.directory) == TRUE) {
    stop("An output directory is needed for the main dataset.")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else {
    dir.create(output.directory)
  }

  if (dir.exists(removed.directory) == TRUE) {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", removed.directory))
      dir.create(removed.directory)
    }
  } else {
    dir.create(removed.directory)
  }
  
  # Gets contig file names
  file.names <- list.files(iupac.directory)

  #Matching and processing for each sample
  for (i in 1:length(file.names)) {

    # Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])

    #Reads in contigs
    contigs = Biostrings::readDNAStringSet(paste0(iupac.directory, "/", file.names[i]), format = "fasta")
    # names(contigs) = paste0("contig_", stringr::str_pad(seq(1:length(contigs)), 6, pad = "0"))

    collect.count = c()
    collect.prop = c()
    for (j in seq_along(contigs)) {
      temp.contig = as.character(contigs[j])
      count = sum(stringr::str_count(temp.contig, pattern = c("R", "Y", "K", "M", "S", "W", "B", "D", "H", "V")))
      temp.prop = count / Biostrings::width(temp.contig)
      collect.count = append(collect.count, count)
      collect.prop = append(collect.count, temp.prop)

    }

    

    # Finds probes that match to two or more contigs
    # final.loci = as.list(as.character(contigs))
    # writeFasta(sequences = final.loci, names = names(final.loci),
    #           paste0(species.dir, "/", sample, "_renamed-contigs.fa"), nbchar = 1000000, as.string = T)

  }#end i loop
    

}#end informative sites function
