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
                            temp.directory = NULL,
                            threads = 1,
                            memory = 1,
                            clean.up = TRUE,
                            overwrite = FALSE,
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

  if (is.null(temp.directory) == TRUE){ temp.directory = getwd() }

  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

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

  mem.cl <- floor(memory / threads)

  # Use mclapply (fork-based) instead of a SOCK cluster to avoid
  # "invalid connection" crashes when a GATK worker process dies.
  parallel::mclapply(seq_along(sample.names), function(i) {

    sample.id      = sample.names[i]
    hap.dir        = paste0(haplotype.caller.directory, "/", sample.id)
    reference.path = paste0(mapping.directory, "/", sample.id, "/index/reference.fa")
    log.file       = paste0("logs/sample_logs/FAILURE_", sample.id, "_baseRecalibration.txt")

    tryCatch({

      gatk = paste0(gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=",
                    temp.directory, " -Xmx", mem.cl, "G\"")

      # Pass 1a: genotype the initial GVCF to get a raw variant set
      system(paste0(gatk, " GenotypeGVCFs -R ", reference.path,
                    " -V ", hap.dir, "/gatk4-haplotype-caller.g.vcf.gz",
                    " --use-new-qual-calculator true",
                    " -O ", hap.dir, "/gatk4-bqsr-genotype.vcf"))

      # Pass 1b: select and hard-filter SNPs
      system(paste0(gatk, " SelectVariants",
                    " -V ", hap.dir, "/gatk4-bqsr-genotype.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-snps.vcf --select-type SNP"))

      system(paste0(gatk, " VariantFiltration -R ", reference.path,
                    " -V ", hap.dir, "/gatk4-bqsr-snps.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-filtered-snps.vcf",
                    " -filter \"QD<2.0\" --filter-name \"QD2\"",
                    " -filter \"QUAL<100.0\" --filter-name \"QUAL30\"",
                    " -filter \"SOR>3.0\" --filter-name \"SOR3\"",
                    " -filter \"FS>60.0\" --filter-name \"FS60\"",
                    " -filter \"MQ<50.0\" --filter-name \"MQ40\"",
                    " -filter \"MQRankSum<-12.5\" --filter-name \"MQRankSum12.5\"",
                    " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum8\""))

      # Pass 1c: select and hard-filter indels
      system(paste0(gatk, " SelectVariants",
                    " -V ", hap.dir, "/gatk4-bqsr-genotype.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-indels.vcf --select-type INDEL"))

      system(paste0(gatk, " VariantFiltration -R ", reference.path,
                    " -V ", hap.dir, "/gatk4-bqsr-indels.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-filtered-indels.vcf",
                    " -filter \"QD<2.0\" --filter-name \"QD2\"",
                    " -filter \"QUAL<100.0\" --filter-name \"QUAL30\"",
                    " -filter \"FS>200.0\" --filter-name \"FS200\"",
                    " -filter \"ReadPosRankSum<-20.0\" --filter-name \"ReadPosRankSum20\""))

      # Pass 1d: merge filtered SNPs + indels into a single "known sites" VCF
      system(paste0(gatk, " SortVcf",
                    " -I ", hap.dir, "/gatk4-bqsr-filtered-snps.vcf",
                    " -I ", hap.dir, "/gatk4-bqsr-filtered-indels.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-filtered-combined.vcf"))

      system(paste0(gatk, " SelectVariants",
                    " -V ", hap.dir, "/gatk4-bqsr-filtered-combined.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-rem-filtered-combined.vcf",
                    " --exclude-filtered TRUE"))

      # Determine BAM path (merged vs. single-lane)
      lane.files = list.dirs(paste0(mapping.directory, "/", sample.id))
      lane.files = lane.files[grep("Lane_", lane.files)]
      read.bam   = if (length(lane.files) == 1) {
        paste0(mapping.directory, "/", sample.id, "/Lane_1")
      } else {
        paste0(mapping.directory, "/", sample.id, "/Lane_Merge")
      }

      # Pass 2a: build recalibration table from the known-sites VCF
      system(paste0(gatk, " BaseRecalibrator",
                    " -I ", read.bam, "/final-mapped-all.bam",
                    " -R ", reference.path,
                    " --known-sites ", hap.dir, "/gatk4-bqsr-rem-filtered-combined.vcf",
                    " -O ", read.bam, "/recal_data.table"))

      # Pass 2b: apply recalibration
      recal.bam = paste0(read.bam, "/bqsr-mapped-all.bam")
      system(paste0(gatk, " ApplyBQSR",
                    " -I ", read.bam, "/final-mapped-all.bam",
                    " -R ", reference.path,
                    " --bqsr-recal-file ", read.bam, "/recal_data.table",
                    " -O ", recal.bam))

      # Pass 2c: re-run HaplotypeCaller on the recalibrated BAM
      system(paste0(gatk, " HaplotypeCaller",
                    " -R ", reference.path,
                    " -I ", recal.bam,
                    " -O ", hap.dir, "/gatk4-bqsr-haplotype-caller.g.vcf.gz",
                    " -ERC GVCF",
                    " -bamout ", hap.dir, "/gatk4-bqsr-haplotype-caller.bam"))

      if (clean.up == TRUE) {
        to.rm = c("gatk4-bqsr-filtered-combined.vcf",
                  "gatk4-bqsr-filtered-indels.vcf",   "gatk4-bqsr-filtered-indels.vcf.idx",
                  "gatk4-bqsr-filtered-snps.vcf",     "gatk4-bqsr-filtered-snps.vcf.idx",
                  "gatk4-bqsr-rem-filtered-combined.vcf", "gatk4-bqsr-rem-filtered-combined.vcf.idx",
                  "gatk4-bqsr-snps.vcf",   "gatk4-bqsr-snps.vcf.idx",
                  "gatk4-bqsr-indels.vcf", "gatk4-bqsr-indels.vcf.idx")
        for (f in to.rm) {
          fp = paste0(hap.dir, "/", f)
          if (file.exists(fp)) { file.remove(fp) }
        }
      }

      print(paste0(sample.id, " completed GATK4 base recalibration!"))

    }, error = function(e) {
      msg = paste0("Unexpected R error: ", conditionMessage(e))
      writeLines(msg, log.file)
      warning(sample.id, ": baseRecalibration failed — see ", log.file)
    })

  }, mc.cores = threads)

  invisible(NULL)

}#end function

# END SCRIPT
