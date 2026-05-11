#' @title filterHeterozygosity
#'
#' @description Filters per-sample FASTA contig files based on the proportion
#'   of IUPAC ambiguity codes (heterozygous sites). Contigs whose proportion of
#'   ambiguous bases equals or exceeds the threshold are moved to a separate
#'   "removed" directory; those below the threshold are written to the output
#'   directory. Samples are processed in parallel.
#'
#' @param iupac.directory path to the directory of input FASTA files containing
#'   IUPAC-coded sequences (e.g. output of VCFtoContigs() with
#'   ambiguity.codes = TRUE).
#'
#' @param output.directory path to the directory where contigs passing the
#'   heterozygosity filter will be saved.
#'
#' @param removed.directory path to the directory where contigs exceeding the
#'   heterozygosity threshold will be saved.
#'
#' @param threshold numeric (0-1); the maximum allowed proportion of IUPAC
#'   ambiguity characters (R, Y, K, M, S, W, B, D, H, V) per contig. Contigs
#'   at or above this proportion are removed.
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total RAM in GB; divided equally across threads.
#'
#' @param overwrite logical; if TRUE existing output and removed directories are
#'   deleted and recreated before processing.
#'
#' @return invisibly; writes filtered FASTA files to output.directory and
#'   high-heterozygosity FASTA files to removed.directory.
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
  cl <- parallel::makeCluster(threads, outfile = "")
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
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

    names(collect.prop) = names(contigs)
    too.high = names(collect.prop[collect.prop >= threshold])

    high.contigs = contigs[names(contigs) %in% too.high]
    low.contigs = contigs[!names(contigs) %in% too.high]

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

  parallel::stopCluster(cl)

} # end function

