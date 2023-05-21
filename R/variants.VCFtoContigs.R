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

VCFtoContigs = function(genotype.directory = "variant-calling",
                        mapping.directory = "sample-mapping",
                        output.directory = "final-contigs",
                        vcf.file = c("SNP", "Indel", "Both"),
                        consensus.sequences = FALSE,
                        ambiguity.codes = TRUE,
                        threads = 1,
                        memory = 1,
                        gatk4.path = NULL,
                        overwrite = TRUE,
                        quiet = TRUE) {

  #Debugging
  # setwd("/Volumes/LaCie/Mantellidae/data-analysis")
  # genotype.directory <- "/Volumes/LaCie/Mantellidae/data-analysis/variant-calling/sample-genotypes"
  # mapping.directory <- "variant-calling/sample-mapping"
  # output.directory = "contigs/5_iupac-contigs"
  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE
  # consensus.sequences = FALSE
  # ambiguity.codes = TRUE
  # vcf.file = "SNP"

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
  if (is.null(genotype.directory) == TRUE) {
    stop("Please provide the genotype directory.")
  }
  if (file.exists(genotype.directory) == FALSE) {
    stop("Genotype directory not found.")
  }

  # Quick checks
  if (is.null(output.directory) == TRUE) {
    stop("Please provide the output directory.")
  }

  if (consensus.sequences == TRUE && ambiguity.codes == TRUE) {
    stop("Cannot have both consensus sequences and ambiguity codes. If both are needed, run function twice")
  }

  vcf.string = NULL
  if (vcf.file == "SNP") {
    vcf.string = "gatk4-final-snps.vcf"
  }

  if (vcf.file == "Indel" || vcf.file == "indel" || vcf.file == "INDEL") {
    vcf.string = "gatk4-final-indels.vcf"
  }

  if (vcf.file == "Both" || vcf.file == "both" || vcf.file == "BOTH") {
    vcf.string <- "gatk4-final-genotypes.vcf"
  }

  if (is.null(vcf.string) == TRUE) {
    stop("please choose SNP, Indel, or Both for vcf.string.")
  }

  # Creates output directory
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  }

  # Get multifile databases together
  sample.names <- list.dirs(genotype.directory, recursive = FALSE, full.names = FALSE)

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
  foreach(i = seq_along(sample.names), .packages = c("foreach", "PhyloProcessR")) %dopar% {

    # Obtains sample vcf
    sample.vcf = paste0(genotype.directory, "/", sample.names[i], "/", vcf.string)
    reference.path = paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")

    if (ambiguity.codes == TRUE) {
      # Selects only the SNPs from the VCF
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
        " FastaAlternateReferenceMaker",
        " -R ", reference.path,
        " -V ", sample.vcf,
        " -O ", output.directory, "/", sample.names[i], ".fa",
        " --use-iupac-sample ", sample.names[i]
      ))
    } # end if
    
    if (consensus.sequences == TRUE) {
      # Selects only the SNPs from the VCF
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
        " FastaAlternateReferenceMaker",
        " -R ", reference.path,
        " -V ", sample.vcf,
        " -O ", output.directory, "/", sample.names[i], ".fa"
      ))
    } # end if

    contigs = Biostrings::readDNAStringSet(paste0(output.directory, "/", sample.names[i], ".fa"), format = "fasta")
    names(contigs) = gsub(":.*", "", names(contigs))
    names(contigs) = gsub(".* ", "", names(contigs))

    system(paste0("rm ", output.directory, "/", sample.names[i], ".fa.fai"))
    system(paste0("rm ", output.directory, "/", sample.names[i], ".dict"))

    # Saves above threshold contigs
    final.loci = as.list(as.character(contigs))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", sample.names[i], ".fa"), nbchar = 1000000, as.string = TRUE
    )

  }#end i loop

  snow::stopCluster(cl)

}#end function

# #END SCRIPT