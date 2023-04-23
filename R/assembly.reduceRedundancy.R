#' @title reduceRedundancy
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param input.reads path to a folder of raw reads in fastq format.
#'
#' @param reference a csv file with a "File" and "Sample" columns, where "File" is the file name and "Sample" is the desired renamed file
#'
#' @param output.name the new directory to save the adaptor trimmed sequences
#'
#' @param mapper "Sample" to run on a single sample or "Directory" to run on a directory of samples
#'
#' @param min.iterations system path to fastp in case it can't be found
#'
#' @param max.iterations system path to fastp in case it can't be found
#'
#' @param min.length system path to fastp in case it can't be found
#'
#' @param max.length system path to fastp in case it can't be found
#'
#' @param min.ref.id system path to fastp in case it can't be found
#'
#' @param spades.path system path to fastp in case it can't be found
#'
#' @param bbmap.path system path to fastp in case it can't be found
#'
#' @param cap3.path system path to fastp in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
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






