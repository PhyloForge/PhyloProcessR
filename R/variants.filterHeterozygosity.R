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
#' astral.data = astralPlane(
#'   astral.tree = your.tree,
#'   outgroups = c("species_one", "species_two"),
#'   tip.length = 1
#' )
#'
#' @export


# Calculates informative sites
filterHeterozygosity = function(iupac.directory = NULL,
                                output.directory = NULL,
                                removed.directory = NULL,
                                threshold = 0.05,
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE) {

  # setwd("/Volumes/LaCie/Mantellidae")
  # iupac.directory = "data-analysis/contigs/5_iupac-contigs"
  # removed.directory = "data-analysis/contigs/6_removed-contigs"
  # output.directory = "data-analysis/contigs/7_filtered-contigs"
  # threshold = 0.05
  # threads = 4
  # memory = 20
  # overwrite = FALSE

  require(foreach)

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


  # Sets up multiprocessing
  cl <- snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl <- floor(memory / threads)

  # Loops through each locus and does operations on them
  foreach(i = seq_along(file.names), .packages = c("foreach", "Biostrings", "stringr", "PhyloProcessR")) %dopar% {

    # Reads in contigs
    contigs = Biostrings::readDNAStringSet(paste0(iupac.directory, "/", file.names[i]), format = "fasta")
    # names(contigs) = paste0("contig_", stringr::str_pad(seq(1:length(contigs)), 6, pad = "0"))

    collect.count = c()
    collect.prop = c()
    for (j in seq_along(contigs)) {
      temp.contig = as.character(contigs[j])
      count = sum(stringr::str_count(temp.contig, pattern = c("R", "Y", "K", "M", "S", "W", "B", "D", "H", "V")))
      temp.prop = count / Biostrings::width(temp.contig)
      collect.count = append(collect.count, count)
      collect.prop = append(collect.prop, temp.prop)
    }

    too.high = collect.prop[collect.prop >= threshold]

    high.contigs = contigs[names(contigs) %in% names(too.high)]
    low.contigs = contigs[!names(contigs) %in% names(too.high)]

    # Saves below threshold contigs
    final.loci = as.list(as.character(low.contigs))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", file.names[i]), nbchar = 1000000, as.string = TRUE
    )

    # Saves above threshold contigs
    final.loci = as.list(as.character(high.contigs))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(removed.directory, "/", file.names[i]), nbchar = 1000000, as.string = TRUE
    )
  } # end i loop

  snow::stopCluster(cl)

} # end function

