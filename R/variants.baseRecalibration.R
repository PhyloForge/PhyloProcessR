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
#' @param temp.directory directory for temporary recalibration files; when
#'   NULL, the current working directory is used.
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total JVM heap budget in GB across concurrent samples. Native
#'   memory and JVM overhead require additional headroom.
#'
#' @param clean.up logical; if TRUE intermediate VCF files generated during the
#'   first-pass genotyping and filtering steps are deleted after the BQSR run.
#'
#' @param overwrite logical; if TRUE samples that already have a BQSR GVCF are
#'   reprocessed.
#'
#' @param quiet logical; currently unused.
#' @param ploidy positive integer ploidy used for the recalibrated caller pass.
#' @param sample.names optional retained sample set.
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
                            quiet = TRUE,
                            ploidy = 2,
                            sample.names = NULL) {

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

  if (length(ploidy) != 1 || !is.finite(ploidy) || ploidy < 1 || ploidy != as.integer(ploidy))
    stop("ploidy must be a positive integer.")

  #Quick checks
  if (is.null(haplotype.caller.directory) == TRUE) {
    stop("Please provide the haplotype caller directory.")
  }
  if (file.exists(haplotype.caller.directory) == F) {
    stop("Haplotype caller directory not found.")
  }

  if (is.null(temp.directory) == TRUE){ temp.directory = tempdir() }
  .ensureDirectory(temp.directory)

  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  #Get multifile databases together
  discovered = list.dirs(haplotype.caller.directory, recursive = F, full.names = F)
  if (is.null(sample.names)) sample.names = discovered

  # Resumes file download
  if (overwrite == FALSE) {
    done.names = sample.names[vapply(sample.names, function(s) {
      d = file.path(haplotype.caller.directory, s)
      g = file.path(d, "gatk4-bqsr-haplotype-caller.g.vcf.gz")
      bam = tryCatch(.selectedSampleBam(mapping.directory, s, TRUE), error = function(e) NA_character_)
      !is.na(bam) && .stageComplete(d, "baseRecalibration", c(g, paste0(g, ".tbi"), bam))
    }, logical(1))]
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0){ return("no samples available to analyze.") }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  resources = .validateResources(threads, memory, length(sample.names))
  gatk.binary = .toolCommand("gatk", gatk4.path)

  # Use mclapply (fork-based) instead of a SOCK cluster to avoid
  # "invalid connection" crashes when a GATK worker process dies.
  results = parallel::mclapply(seq_along(sample.names), function(i) {

    sample.id      = sample.names[i]
    hap.dir        = paste0(haplotype.caller.directory, "/", sample.id)
    reference.path = paste0(mapping.directory, "/", sample.id, "/index/reference.fa")
    log.file       = paste0("logs/sample_logs/FAILURE_", sample.id, "_baseRecalibration.txt")

    tryCatch({

      gatk = .gatkCommand(gatk.binary, temp.directory, resources$heap.mb)

      # Pass 1a: genotype the initial GVCF to get a raw variant set
      .runCommand(paste0(gatk, " GenotypeGVCFs -R ", shQuote(reference.path),
                    " -V ", hap.dir, "/gatk4-haplotype-caller.g.vcf.gz",
                    " --use-new-qual-calculator true",
                    " -O ", shQuote(file.path(hap.dir, "gatk4-bqsr-genotype.vcf"))), quiet, "BQSR genotyping", stderr.log = log.file)

      # Pass 1b: select and hard-filter SNPs
      .runCommand(paste0(gatk, " SelectVariants",
                    " -V ", hap.dir, "/gatk4-bqsr-genotype.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-snps.vcf --select-type SNP"), quiet, "BQSR SNP selection", stderr.log = log.file)

      .runCommand(paste0(gatk, " VariantFiltration -R ", reference.path,
                    " -V ", hap.dir, "/gatk4-bqsr-snps.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-filtered-snps.vcf",
                    " -filter \"QD<2.0\" --filter-name \"QD2\"",
                    " -filter \"QUAL<100.0\" --filter-name \"QUAL100\"",
                    " -filter \"SOR>3.0\" --filter-name \"SOR3\"",
                    " -filter \"FS>60.0\" --filter-name \"FS60\"",
                    " -filter \"MQ<50.0\" --filter-name \"MQ50\"",
                    " -filter \"MQRankSum<-12.5\" --filter-name \"MQRankSum12.5\"",
                    " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum8\""), quiet, "BQSR SNP filtering", stderr.log = log.file)

      # Pass 1c: select and hard-filter indels
      .runCommand(paste0(gatk, " SelectVariants",
                    " -V ", hap.dir, "/gatk4-bqsr-genotype.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-indels.vcf --select-type INDEL"), quiet, "BQSR indel selection", stderr.log = log.file)

      .runCommand(paste0(gatk, " VariantFiltration -R ", reference.path,
                    " -V ", hap.dir, "/gatk4-bqsr-indels.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-filtered-indels.vcf",
                    " -filter \"QD<2.0\" --filter-name \"QD2\"",
                    " -filter \"QUAL<100.0\" --filter-name \"QUAL100\"",
                    " -filter \"FS>200.0\" --filter-name \"FS200\"",
                    " -filter \"ReadPosRankSum<-20.0\" --filter-name \"ReadPosRankSum20\""), quiet, "BQSR indel filtering", stderr.log = log.file)

      # Pass 1d: merge filtered SNPs + indels into a single "known sites" VCF
      .runCommand(paste0(gatk, " SortVcf",
                    " -I ", hap.dir, "/gatk4-bqsr-filtered-snps.vcf",
                    " -I ", hap.dir, "/gatk4-bqsr-filtered-indels.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-filtered-combined.vcf"), quiet, "BQSR VCF merge", stderr.log = log.file)

      .runCommand(paste0(gatk, " SelectVariants",
                    " -V ", hap.dir, "/gatk4-bqsr-filtered-combined.vcf",
                    " -O ", hap.dir, "/gatk4-bqsr-rem-filtered-combined.vcf",
                    " --exclude-filtered TRUE"), quiet, "BQSR passing variant selection", stderr.log = log.file)

      # Determine BAM path (merged vs. single-lane)
      source.bam = .selectedSampleBam(mapping.directory, sample.id, FALSE)
      read.bam = dirname(source.bam)

      # Pass 2a: build recalibration table from the known-sites VCF
      .runCommand(paste0(gatk, " BaseRecalibrator",
                    " -I ", source.bam,
                    " -R ", reference.path,
                    " --known-sites ", hap.dir, "/gatk4-bqsr-rem-filtered-combined.vcf",
                    " -O ", read.bam, "/recal_data.table"), quiet, "BaseRecalibrator", stderr.log = log.file)

      # Pass 2b: apply recalibration
      recal.bam = paste0(read.bam, "/bqsr-mapped-all.bam")
      .runCommand(paste0(gatk, " ApplyBQSR",
                    " -I ", source.bam,
                    " -R ", reference.path,
                    " --bqsr-recal-file ", read.bam, "/recal_data.table",
                    " -O ", recal.bam), quiet, "ApplyBQSR", stderr.log = log.file)

      # Pass 2c: re-run HaplotypeCaller on the recalibrated BAM
      .runCommand(paste0(gatk, " HaplotypeCaller",
                    " -R ", reference.path,
                    " -I ", recal.bam,
                    " -O ", hap.dir, "/gatk4-bqsr-haplotype-caller.g.vcf.gz",
                    " -ERC GVCF -ploidy ", ploidy,
                    " --native-pair-hmm-threads 1 -bamout ", hap.dir, "/gatk4-bqsr-haplotype-caller.bam"), quiet, "BQSR HaplotypeCaller", stderr.log = log.file)

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

      output = file.path(hap.dir, "gatk4-bqsr-haplotype-caller.g.vcf.gz")
      if (!all(file.exists(c(output, paste0(output, ".tbi"), recal.bam)))) stop("BQSR outputs are incomplete")
      .markStageComplete(hap.dir, "baseRecalibration", c(paste0("ploidy=", ploidy), paste0("bam=", source.bam)))
      list(success = TRUE)

    }, error = function(e) {
      msg = paste0("Unexpected R error: ", conditionMessage(e))
      cat("\n", msg, "\n", file = log.file, append = TRUE)
      list(success = FALSE, message = conditionMessage(e))
    })

  }, mc.cores = resources$workers)

  .collectWorkers(results, sample.names, "Base recalibration")
  invisible(sample.names)

}#end function

# END SCRIPT
