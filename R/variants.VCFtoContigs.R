#' @title VCFtoContigs
#'
#' @description Converts per-sample VCF variant files into FASTA contig
#'   sequences using GATK4 FastaAlternateReferenceMaker. For each sample the
#'   chosen VCF (SNPs, indels, or both) is applied to the sample's reference
#'   FASTA to produce either IUPAC ambiguity-coded sequences or straight
#'   consensus sequences. Samples are processed in parallel.
#'
#' @param genotype.directory path to the directory containing per-sample
#'   genotype sub-directories (output of genotypeSamples()).
#'
#' @param mapping.directory path to the directory containing per-sample
#'   reference index sub-directories (output of mapReferenceSample() or
#'   mapReferenceConsensus()), used to locate each sample's reference.fa.
#'
#' @param output.directory path to the directory where output FASTA contig
#'   files (one per sample) will be saved.
#'
#' @param vcf.file character; which variant type to apply: "SNP" uses the
#'   SNP-only VCF, "Indel" uses the indel-only VCF, "Both" uses the combined
#'   genotype VCF.
#'
#' @param consensus.sequences logical; if TRUE GATK produces a consensus
#'   sequence (majority allele). Cannot be TRUE when ambiguity.codes is TRUE.
#'
#' @param ambiguity.codes logical; if TRUE GATK produces IUPAC ambiguity codes
#'   at heterozygous sites (--use-iupac-sample). Cannot be TRUE when
#'   consensus.sequences is TRUE.
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total RAM in GB; divided equally across threads.
#'
#' @param temp.directory path to a dedicated directory for GATK JVM and Picard
#'   sorting temp files. NULL defaults to the current working directory, which
#'   will cause temp files to accumulate in the project root.
#'
#' @param gatk4.path system path to the directory containing the gatk
#'   executable; NULL searches the system PATH.
#'
#' @param overwrite logical; if TRUE already-completed samples are reprocessed.
#'
#' @param quiet logical; currently unused.
#'
#' @return invisibly; writes one FASTA file per sample to output.directory.
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
                        temp.directory = NULL,
                        gatk4.path = NULL,
                        overwrite = FALSE,
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

  if (is.null(temp.directory) == TRUE) { temp.directory = getwd() }

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
    done.files <- list.files(output.directory, full.names = FALSE)
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

  mem.cl <- floor(memory / threads)

  #Loops through each locus and does operations on them
  parallel::mclapply(seq_along(sample.names), function(i) {
  tryCatch({

    # Obtains sample vcf
    sample.vcf = paste0(genotype.directory, "/", sample.names[i], "/", vcf.string)
    reference.path = paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")

    gatk = paste0(gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=",
                  temp.directory, " -Xmx", mem.cl, "G\"")

    if (ambiguity.codes == TRUE) {
      system(paste0(gatk, " FastaAlternateReferenceMaker",
                    " -R ", reference.path,
                    " -V ", sample.vcf,
                    " -O ", output.directory, "/", sample.names[i], ".fa",
                    " --use-iupac-sample ", sample.names[i]))
    }

    if (consensus.sequences == TRUE) {
      system(paste0(gatk, " FastaAlternateReferenceMaker",
                    " -R ", reference.path,
                    " -V ", sample.vcf,
                    " -O ", output.directory, "/", sample.names[i], ".fa"))
    }

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

  }, error = function(e) {
    warning(sample.names[i], " failed: ", conditionMessage(e))
  })
  }, mc.cores = threads) #end i loop

}#end function

# #END SCRIPT