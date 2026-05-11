#' @title reduceRedundancy
#'
#' @description Clusters and de-replicates contigs within each sample assembly
#'   using \code{cd-hit-est} at a user-specified sequence identity threshold.
#'   Each FASTA file in the input directory is processed in parallel: cd-hit-est
#'   collapses highly similar sequences, contigs are renamed sequentially, and
#'   the dereplicated assembly is written to the output directory. The cluster
#'   report (\code{.clstr}) files are removed after processing.
#'
#' @param assembly.directory path to a directory of assembly FASTA files, one
#'   per sample.
#'
#' @param output.directory path to the directory where dereplicated FASTA files
#'   will be written. Default: \code{"reduced-redundancy"}.
#'
#' @param similarity sequence identity threshold for cd-hit-est clustering
#'   (0.6–1.0). Must be >= 0.6. Default: \code{0.95}.
#'
#' @param cdhit.path path to the directory containing \code{cd-hit-est}. If
#'   \code{NULL} expected on the system PATH. Default: \code{NULL}.
#'
#' @param memory total RAM in GB to divide across parallel workers. Default:
#'   \code{1}.
#'
#' @param threads number of parallel workers. Default: \code{1}.
#'
#' @param overwrite logical; if \code{TRUE} the output directory is deleted and
#'   recreated. Default: \code{FALSE}.
#'
#' @param quiet logical; reserved for future use. Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Writes one dereplicated FASTA file per
#'   sample to \code{output.directory}, with contigs renamed
#'   \code{contig_1}, \code{contig_2}, etc.
#'
#' @export

#Iteratively assembles to reference
reduceRedundancy = function(assembly.directory = NULL,
                            output.directory = "reduced-redundancy",
                            similarity = 0.95,
                            cdhit.path = NULL,
                            memory = 1,
                            threads = 1,
                            overwrite = FALSE,
                            quiet = TRUE) {

  # #Debug
  # library(PhyloProcessR)
  # setwd("/Volumes/LaCie/Mantellidae")
  # assembly.directory = "data-analysis/contigs/draft-assemblies"
  # output.directory = "data-analysis/contigs/reduced-redundancy"
  # cdhit.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # similarity <- 0.95
  # quiet = TRUE
  # overwrite = TRUE
  # threads = 5
  # memory = 30

  require(foreach)

  # Same adds to bbmap path
  if (is.null(cdhit.path) == FALSE) {
    b.string <- unlist(strsplit(cdhit.path, ""))
    if (b.string[length(b.string)] != "/") {
      cdhit.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    cdhit.path <- ""
  }

  # Quick checks
  if (is.null(assembly.directory) == TRUE) {
    stop("Please provide input reads.")
  }
  if (is.null(output.directory) == TRUE) {
    stop("Please provide an output directory.")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else {
    dir.create(output.directory)
  }

  file.names = list.files(assembly.directory)
  
  if (similarity >= 0.7) {
    n.val = 5
  }

  if (similarity < 0.7 && similarity >= 0.6) {
    n.val <- 4
  }

  if (similarity < 0.6) {
    stop("similarity too small for cd-hit est. Must be greater than 0.6")
  }
  
  # Sets up multiprocessing
  cl <- snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl <- floor(memory / threads)
  
  #Loop for cd-hit est reductions
  foreach::foreach(i = seq_along(file.names), .packages = c("foreach", "PhyloProcessR")) %dopar% {

    system(paste0(
      cdhit.path, "cd-hit-est -i ", assembly.directory, "/", file.names[i],
      " -o ", output.directory, "/", file.names[i], " -p 0 -T 1",
      " -n ", n.val, " -c ", similarity, " -M ", mem.cl * 100
    ))

    ### Read in data
    all.data = Biostrings::readDNAStringSet(file = paste0(output.directory, "/", file.names[i]), format = "fasta")

    names(all.data) = paste0("contig_", seq(seq_along(all.data)))

    # Writes the final loci
    final.loci = as.list(as.character(all.data))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", file.names[i]),
      nbchar = 1000000, as.string = TRUE, open = "w"
    )

    system(paste0(" rm ", output.directory, "/", file.names[i], ".clstr"))

  }#end iterations if

snow::stopCluster(cl)

  ##########################
}#end function






