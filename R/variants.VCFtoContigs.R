#' @title VCFtoContigs
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param read.directory directory of processed reads
#'
#' @param output.directory save name for the output directory
#'
#' @param full.path.spades contigs are added into existing alignment if algorithm is "add"
#'
#' @param mismatch.corrector algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param kmer.values if a file name is provided, save.name will be used to save aligment to file as a fasta
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

#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
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

VCFtoContigs = function(haplotype.caller.directory = "variant-calling/haplotype-caller",
                        sample.mapping.directory = "variant-calling/sample-mapping",
                        output.directory = "variant-calling/final-contigs",
                        invariant.sites = TRUE,
                        ambiguity.codes = FALSE,
                        haplotype.split = FALSE,
                        save.all.types = FALSE,
                        threads = 1,
                        memory = 1,
                        overwrite = TRUE,
                        quiet = TRUE) {

  # #Debugging
  # haplotype.caller.directory <- "/Volumes/LaCie/Mantellidae/data-analysis/variant-calling/haplotype-caller"
  # sample.mapping.directory <- "variant-calling/sample-mapping"
  # output.directory = "variant-calling/final-contigs"
  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE

  # Same adds to bbmap path
  require(foreach)

  # Same adds to bbmap path
  if (is.null(gatk4.path) == FALSE) {
    b.string <- unlist(strsplit(gatk4.path, ""))
    if (b.string[length(b.string)] != "/") {
      gatk4.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    gatk4.path <- ""
  }

  # Quick checks
  if (is.null(haplotype.caller.directory) == TRUE) {
    stop("Please provide the haplotype caller directory.")
  }
  if (file.exists(haplotype.caller.directory) == FALSE) {
    stop("Haplotype caller directory not found.")
  }

  # Quick checks
  if (is.null(sample.mapping.directory) == TRUE) {
    stop("Please provide the sample mapping directory.")
  }
  if (file.exists(sample.mapping.directory) == FALSE) {
    stop("Sample mapping directory not found.")
  }


  # Creates output directory
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  }

  # Get multifile databases together
  sample.names <- list.dirs(haplotype.caller.directory, recursive = FALSE, full.names = FALSE)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(output.directory, full.names = TRUE, recursive = TRUE)
    done.names <- gsub(".fa$", "", done.files)
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0) {
    return("no samples available to analyze.")
  }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  #Sets up multiprocessing
  cl <- snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl <- floor(memory / threads)

  #Loops through each locus and does operations on them
  foreach(i = seq_along(sample.names), .packages = c("foreach", "vcfR", "PhyloProcessR", "ape", "Biostrings")) %dopar% {

    # Obtains sample vcf
    sample.vcf = paste0(haplotype.caller.directory, "/", sample.names[i], "/gatk4-final-genotype.vcf")
    reference.path = paste0(sample.mapping.directory, "/", sample.names[i], "/index/reference.fa")

    # Selects only the SNPs from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " FastaAlternateReferenceMaker",
      " -R ", reference.path,
      " -V ", sample.vcf,
      " -O ", output.directory, "/", sample.names[i], ".fa",
      " --use-iupac-sample ", sample.names[i]
    ))

  }#end i loop

  snow::stopCluster(cl)

}#end function

# #END SCRIPT