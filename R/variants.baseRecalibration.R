#' @title baseRecalibration
#'
#' @description Performs GATK4 base quality score recalibration (BQSR) as a
#'   two-pass procedure. In pass 1, an initial round of genotyping and hard
#'   filtering is run on each sample's existing haplotype caller GVCF to
#'   produce a set of "known variants". In pass 2, GATK BaseRecalibrator and
#'   ApplyBQSR use those variants to recalibrate the base quality scores in the
#'   original BAM file. A final HaplotypeCaller run then produces a new GVCF
#'   from the recalibrated BAM. Samples are processed in parallel.
#'
#' @param haplotype.caller.directory path to the directory containing per-sample
#'   sub-directories of initial haplotype caller GVCF files (output of
#'   haplotypeCaller()).
#'
#' @param mapping.directory path to the directory containing per-sample BAM
#'   files and reference indices (output of mapReferenceSample() or
#'   mapReferenceConsensus()).
#'
#' @param gatk4.path system path to the directory containing the gatk
#'   executable; NULL searches the system PATH.
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total RAM in GB to allocate; used per-thread as a JVM heap
#'   (-Xmx).
#'
#' @param clean.up logical; if TRUE intermediate VCF files generated during the
#'   first-pass genotyping and filtering steps are deleted after the BQSR run.
#'
#' @param overwrite logical; if TRUE samples that already have a BQSR GVCF are
#'   reprocessed.
#'
#' @param quiet logical; currently unused.
#'
#' @return invisibly; writes recalibrated GVCFs to haplotype.caller.directory
#'   and recalibrated BAMs to mapping.directory.
#'
#' @export

baseRecalibration = function(haplotype.caller.directory = "haplotype-caller",
                            mapping.directory = "sample-mapping",
                            gatk4.path = NULL,
                            threads = 1,
                            memory = 1,
                            clean.up = TRUE,
                            overwrite = TRUE,
                            quiet = TRUE) {

  #Debugging
  #Home directoroies
  # library(PhyloCap)
  # setwd("/Volumes/LaCie/Mantellidae/data-analysis")
  # haplotype.caller.directory <- "variant-calling/haplotype-caller"
  # mapping.directory <- "variant-calling/sample-mapping"

  # gatk4.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"

  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE
  # clean.up = TRUE

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

  #Quick checks
  if (is.null(haplotype.caller.directory) == TRUE) {
    stop("Please provide the haplotype caller directory.")
  }
  if (file.exists(haplotype.caller.directory) == F) {
    stop("Haplotype caller directory not found.")
  }

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Get multifile databases together
  sample.names = list.dirs(haplotype.caller.directory, recursive = F, full.names = F)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(haplotype.caller.directory, full.names = T, recursive = T)
    done.files <- done.files[grep("gatk4-bqsr-haplotype-caller.g.vcf.gz$", done.files)]
    done.names <- gsub("/gatk4-bqsr-haplotype-caller.g.vcf.gz$", "", done.files)
    done.names <- gsub(".*\\/", "", done.names)
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0){ return("no samples available to analyze.") }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  # Sets up multiprocessing
  cl <- snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl <- floor(memory / threads)

  #Loops through each locus and does operations on them
  foreach::foreach(i=seq_along(sample.names), .packages = c("foreach")) %dopar% {
    # Loops through each locus and does operations on them
    # for (i in 1:length(loci.names)){

    reference.path <- paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " GenotypeGVCFs -R ", reference.path,
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-haplotype-caller.g.vcf.gz",
      " --use-new-qual-calculator true -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-genotype.vcf"
    ))

    # Selects only the SNPs from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-genotype.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf",
      " --select-type SNP"
    ))

    #Selects only the indels from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-genotype.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf",
      " --select-type INDEL"
    ))

    # Applies filters to SNPs
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf",
      " -filter \"QD<2.0\" --filter-name \"QD2\"",
      " -filter \"QUAL<100.0\" --filter-name \"QUAL30\"",
      " -filter \"SOR>3.0\" --filter-name \"SOR3\"",
      " -filter \"FS>60.0\" --filter-name \"FS60\"",
      " -filter \"MQ<50.0\" --filter-name \"MQ40\"",
      " -filter \"MQRankSum<-12.5\" --filter-name \"MQRankSum12.5\"",
      " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum8\""
    ))

    # Applies filters to indels
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf",
      " -filter \"QD<2.0\" --filter-name \"QD2\"",
      " -filter \"QUAL<100.0\" --filter-name \"QUAL30\"",
      " -filter \"FS>200.0\" --filter-name \"FS200\"",
      " -filter \"ReadPosRankSum<-20.0\" --filter-name \"ReadPosRankSum20\""
    ))

    # Combine them into a single VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SortVcf",
      " -I ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf",
      " -I ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-combined.vcf"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-combined.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf",
      " --exclude-filtered TRUE"
    ))

    #Gathers correct read files
    lane.files = list.dirs(paste0(mapping.directory, "/", sample.names[i]))
    lane.files = lane.files[grep("Lane_", lane.files)]

    if (length(lane.files) != 1) {
      read.bam = paste0(mapping.directory, "/", sample.names[i], "/Lane_Merge")
    }
    if (length(lane.files) == 1) {
      read.bam = paste0(mapping.directory, "/", sample.names[i], "/Lane_1")
    }

    #Applies base recalibration to the original bam file
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " BaseRecalibrator",
      " -I ", read.bam, "/final-mapped-all.bam",
      " -R ", reference.path,
      " --known-sites ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf",
      " -O ", read.bam, "/recal_data.table"
    ))

    recal.bam <- paste0(read.bam, "/bqsr-mapped-all.bam")

    # Applies the recalibration and creates a new bam file
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " ApplyBQSR",
      " -I ", read.bam, "/final-mapped-all.bam",
      " -R ", reference.path,
      " --bqsr-recal-file ", read.bam, "/recal_data.table",
      " -O ", recal.bam
    ))

   # Starts to finally look for Haplotypes! *here
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " HaplotypeCaller",
      " -R ", reference.path,
      " -I ", recal.bam,
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-haplotype-caller.g.vcf.gz",
      " -ERC GVCF",
      " -bamout ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-haplotype-caller.bam"
    ))

    if (clean.up == TRUE){

      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-combined.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf.idx"))

    }

    print(paste0(sample.names[i], " completed GATK4 base recalibration!"))

  }#end i loop

  snow::stopCluster(cl)

}#end function

# END SCRIPT
