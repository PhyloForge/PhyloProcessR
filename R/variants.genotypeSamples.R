#' Genotype and hard-filter per-sample GVCFs
#' @param mapping.directory Mapping directory.
#' @param haplotype.caller.directory GVCF directory.
#' @param output.directory Genotype output directory.
#' @param use.base.recalibration Select BQSR GVCFs.
#' @param custom.SNP.QD,custom.SNP.QUAL,custom.SNP.SOR,custom.SNP.FS,custom.SNP.MQ,custom.SNP.MQRankSum,custom.SNP.ReadPosRankSum SNP record thresholds.
#' @param custom.INDEL.QD,custom.INDEL.QUAL,custom.INDEL.FS,custom.INDEL.ReadPosRankSum indel record thresholds.
#' @param gatk4.path GATK path/directory or NULL.
#' @param temp.directory Temporary directory.
#' @param threads Concurrent samples.
#' @param memory Total JVM heap budget in GB.
#' @param overwrite Recompute outputs.
#' @param quiet Suppress tool output while retaining logs.
#' @param sample.names Optional retained sample set.
#' @return Invisibly returns sample names.
#' @export
genotypeSamples = function(mapping.directory = NULL,
                           haplotype.caller.directory = NULL,
                           output.directory = "sample-genotypes",
                           use.base.recalibration = FALSE,
                           custom.SNP.QD = 2,
                           custom.SNP.QUAL = 30,
                           custom.SNP.SOR = 3,
                           custom.SNP.FS = 60,
                           custom.SNP.MQ = 40,
                           custom.SNP.MQRankSum = -12.5,
                           custom.SNP.ReadPosRankSum = -8,
                           custom.INDEL.QD = 2,
                           custom.INDEL.QUAL = 30,
                           custom.INDEL.FS = 60,
                           custom.INDEL.ReadPosRankSum = -8,
                           gatk4.path = NULL,
                           temp.directory = NULL,
                           threads = 1,
                           memory = 1,
                           overwrite = FALSE,
                           quiet = TRUE,
                           sample.names = NULL) {

  #Quick checks: every threshold must be one finite number
  thresholds = list(custom.SNP.QD = custom.SNP.QD,
                    custom.SNP.QUAL = custom.SNP.QUAL,
                    custom.SNP.SOR = custom.SNP.SOR,
                    custom.SNP.FS = custom.SNP.FS,
                    custom.SNP.MQ = custom.SNP.MQ,
                    custom.SNP.MQRankSum = custom.SNP.MQRankSum,
                    custom.SNP.ReadPosRankSum = custom.SNP.ReadPosRankSum,
                    custom.INDEL.QD = custom.INDEL.QD,
                    custom.INDEL.QUAL = custom.INDEL.QUAL,
                    custom.INDEL.FS = custom.INDEL.FS,
                    custom.INDEL.ReadPosRankSum = custom.INDEL.ReadPosRankSum)
  bad = names(thresholds)[!vapply(thresholds, function(x) {
    length(x) == 1 && is.numeric(x) && is.finite(x)
  }, logical(1))]
  if (length(bad) > 0) { stop(bad[1], " must be one finite numeric value.") }

  if (is.null(mapping.directory) || !dir.exists(mapping.directory)) {
    stop("Mapping directory not found.")
  }
  if (is.null(haplotype.caller.directory) || !dir.exists(haplotype.caller.directory)) {
    stop("Haplotype caller directory not found.")
  }

  discovered = list.dirs(haplotype.caller.directory, recursive = FALSE, full.names = FALSE)
  if (is.null(sample.names)) { sample.names = discovered }
  if (length(sample.names) == 0) { stop("No samples are available to genotype.") }

  #Selects the initial or the recalibrated GVCF, and the per-sample reference
  gvcf.name = if (use.base.recalibration == TRUE) {
    "gatk4-bqsr-haplotype-caller.g.vcf.gz"
  } else {
    "gatk4-haplotype-caller.g.vcf.gz"
  }
  inputs = file.path(haplotype.caller.directory, sample.names, gvcf.name)
  refs = file.path(mapping.directory, sample.names, "index", "reference.fa")
  if (any(!file.exists(inputs))) {
    stop("Missing selected GVCF for sample(s): ",
         paste(sample.names[!file.exists(inputs)], collapse = ", "))
  }
  if (any(!file.exists(refs))) {
    stop("Missing reference for sample(s): ",
         paste(sample.names[!file.exists(refs)], collapse = ", "))
  }

  resources = .validateResources(threads, memory, length(sample.names))
  gatk = .toolCommand("gatk", gatk4.path)
  if (is.null(temp.directory)) { temp.directory = tempdir() }
  .ensureDirectory(temp.directory)
  if (overwrite == TRUE && dir.exists(output.directory)) {
    unlink(output.directory, recursive = TRUE)
  }
  .ensureDirectory(output.directory)
  .ensureDirectory("logs/sample_logs")
  command = .gatkCommand(gatk, temp.directory, resources$heap.mb)

  #Genotypes each sample in parallel
  results = parallel::mclapply(seq_along(sample.names), function(i) {
    tryCatch({
      sample = sample.names[i]
      sample.dir = file.path(output.directory, sample)
      .ensureDirectory(sample.dir)

      final = file.path(sample.dir, c("gatk4-final-genotypes.vcf",
                                      "gatk4-final-snps.vcf",
                                      "gatk4-final-indels.vcf"))
      indices = paste0(final, ".idx")
      if (overwrite == FALSE && .stageComplete(sample.dir, "genotypeSamples", c(final, indices))) {
        return(list(success = TRUE))
      }
      .invalidateStage(sample.dir, "genotypeSamples")
      file.remove(c(final, indices))
      log = file.path("logs/sample_logs", sample, "genotypeSamples.stderr.log")
      .ensureDirectory(dirname(log))

      #Intermediate VCF paths
      unfiltered = file.path(sample.dir, "gatk4-unfiltered-genotypes.vcf")
      unfiltered.snps = file.path(sample.dir, "gatk4-unfiltered-snps.vcf")
      unfiltered.indels = file.path(sample.dir, "gatk4-unfiltered-indels.vcf")
      filtered.snps = file.path(sample.dir, "gatk4-filtered-snps.vcf")
      filtered.indels = file.path(sample.dir, "gatk4-filtered-indels.vcf")
      filtered.genotypes = file.path(sample.dir, "gatk4-filtered-genotypes.vcf")

      #Hard-filter expressions: a record failing any test is labeled and removed
      snp.expr = paste0("QD < ", custom.SNP.QD,
                        " || QUAL < ", custom.SNP.QUAL,
                        " || SOR > ", custom.SNP.SOR,
                        " || FS > ", custom.SNP.FS,
                        " || MQ < ", custom.SNP.MQ,
                        " || MQRankSum < ", custom.SNP.MQRankSum,
                        " || ReadPosRankSum < ", custom.SNP.ReadPosRankSum)
      indel.expr = paste0("QD < ", custom.INDEL.QD,
                          " || QUAL < ", custom.INDEL.QUAL,
                          " || FS > ", custom.INDEL.FS,
                          " || ReadPosRankSum < ", custom.INDEL.ReadPosRankSum)

      #The genotyping steps run in order: genotype, split SNPs and indels, hard
      #filter each, merge, then keep only the passing records.
      cmds = c(
        paste(command, "GenotypeGVCFs -R", shQuote(refs[i]), "-V", shQuote(inputs[i]),
              "--use-new-qual-calculator true -O", shQuote(unfiltered)),
        paste(command, "SelectVariants -R", shQuote(refs[i]), "-V", shQuote(unfiltered),
              "-O", shQuote(unfiltered.snps), "--select-type-to-include SNP"),
        paste(command, "SelectVariants -R", shQuote(refs[i]), "-V", shQuote(unfiltered),
              "-O", shQuote(unfiltered.indels), "--select-type-to-include INDEL"),
        paste(command, "VariantFiltration -R", shQuote(refs[i]), "-V", shQuote(unfiltered.snps),
              "-O", shQuote(filtered.snps),
              "--filter-name SNP_HARD_FILTER --filter-expression", shQuote(snp.expr)),
        paste(command, "VariantFiltration -R", shQuote(refs[i]), "-V", shQuote(unfiltered.indels),
              "-O", shQuote(filtered.indels),
              "--filter-name INDEL_HARD_FILTER --filter-expression", shQuote(indel.expr)),
        paste(command, "MergeVcfs -I", shQuote(filtered.snps), "-I", shQuote(filtered.indels),
              "-O", shQuote(filtered.genotypes)),
        paste(command, "SelectVariants -R", shQuote(refs[i]), "-V", shQuote(filtered.genotypes),
              "-O", shQuote(final[1]), "--exclude-filtered true"),
        paste(command, "SelectVariants -R", shQuote(refs[i]), "-V", shQuote(filtered.snps),
              "-O", shQuote(final[2]), "--exclude-filtered true"),
        paste(command, "SelectVariants -R", shQuote(refs[i]), "-V", shQuote(filtered.indels),
              "-O", shQuote(final[3]), "--exclude-filtered true"))
      for (k in seq_along(cmds)) {
        .runCommand(cmds[k], quiet, paste0("genotyping command ", k), stderr.log = log)
      }

      if (!all(file.exists(c(final, indices)))) {
        stop("Final VCF set or indices were not created")
      }
      .markStageComplete(sample.dir, "genotypeSamples")
      list(success = TRUE)
    }, error = function(e) list(success = FALSE, message = conditionMessage(e)))
  }, mc.cores = resources$workers)

  .collectWorkers(results, sample.names, "Genotyping")
  invisible(sample.names)
}#end function
