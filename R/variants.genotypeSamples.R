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
genotypeSamples = function(mapping.directory = NULL, haplotype.caller.directory = NULL,
 output.directory = "sample-genotypes", use.base.recalibration = FALSE,
 custom.SNP.QD = 2, custom.SNP.QUAL = 30, custom.SNP.SOR = 3, custom.SNP.FS = 60,
 custom.SNP.MQ = 40, custom.SNP.MQRankSum = -12.5, custom.SNP.ReadPosRankSum = -8,
 custom.INDEL.QD = 2, custom.INDEL.QUAL = 30, custom.INDEL.FS = 60,
 custom.INDEL.ReadPosRankSum = -8, gatk4.path = NULL, temp.directory = NULL,
 threads = 1, memory = 1, overwrite = FALSE, quiet = TRUE, sample.names = NULL) {
  thresholds = list(custom.SNP.QD=custom.SNP.QD, custom.SNP.QUAL=custom.SNP.QUAL,
    custom.SNP.SOR=custom.SNP.SOR, custom.SNP.FS=custom.SNP.FS, custom.SNP.MQ=custom.SNP.MQ,
    custom.SNP.MQRankSum=custom.SNP.MQRankSum, custom.SNP.ReadPosRankSum=custom.SNP.ReadPosRankSum,
    custom.INDEL.QD=custom.INDEL.QD, custom.INDEL.QUAL=custom.INDEL.QUAL,
    custom.INDEL.FS=custom.INDEL.FS, custom.INDEL.ReadPosRankSum=custom.INDEL.ReadPosRankSum)
  bad = names(thresholds)[!vapply(thresholds, function(x) length(x)==1 && is.numeric(x) && is.finite(x), logical(1))]
  if (length(bad)) stop(bad[1], " must be one finite numeric value.")
  if (is.null(mapping.directory) || !dir.exists(mapping.directory)) stop("Mapping directory not found.")
  if (is.null(haplotype.caller.directory) || !dir.exists(haplotype.caller.directory)) stop("Haplotype caller directory not found.")
  discovered = list.dirs(haplotype.caller.directory, recursive = FALSE, full.names = FALSE)
  if (is.null(sample.names)) sample.names = discovered
  if (!length(sample.names)) stop("No samples are available to genotype.")
  gname = if (use.base.recalibration) "gatk4-bqsr-haplotype-caller.g.vcf.gz" else "gatk4-haplotype-caller.g.vcf.gz"
  inputs = file.path(haplotype.caller.directory, sample.names, gname)
  refs = file.path(mapping.directory, sample.names, "index", "reference.fa")
  if (any(!file.exists(inputs))) stop("Missing selected GVCF for sample(s): ", paste(sample.names[!file.exists(inputs)], collapse = ", "))
  if (any(!file.exists(refs))) stop("Missing reference for sample(s): ", paste(sample.names[!file.exists(refs)], collapse = ", "))
  resources = .validateResources(threads, memory, length(sample.names)); gatk = .toolCommand("gatk", gatk4.path)
  if (is.null(temp.directory)) temp.directory = tempdir(); .ensureDirectory(temp.directory)
  if (overwrite && dir.exists(output.directory)) unlink(output.directory, recursive = TRUE)
  .ensureDirectory(output.directory); .ensureDirectory("logs/sample_logs")
  command = .gatkCommand(gatk, temp.directory, resources$heap.mb)
  results = parallel::mclapply(seq_along(sample.names), function(i) tryCatch({
    s=sample.names[i]; d=file.path(output.directory,s); .ensureDirectory(d)
    final = file.path(d,c("gatk4-final-genotypes.vcf","gatk4-final-snps.vcf","gatk4-final-indels.vcf")); indices=paste0(final,".idx")
    if (!overwrite && .stageComplete(d,"genotypeSamples",c(final,indices))) return(list(success=TRUE))
    .invalidateStage(d,"genotypeSamples"); file.remove(c(final,indices)); log=file.path("logs/sample_logs",s,"genotypeSamples.stderr.log"); .ensureDirectory(dirname(log))
    f=function(n) file.path(d,n); unfiltered=f("gatk4-unfiltered-genotypes.vcf"); us=f("gatk4-unfiltered-snps.vcf"); ui=f("gatk4-unfiltered-indels.vcf"); fs=f("gatk4-filtered-snps.vcf"); fi=f("gatk4-filtered-indels.vcf"); fg=f("gatk4-filtered-genotypes.vcf")
    snp.expr=paste0("QD < ",custom.SNP.QD," || QUAL < ",custom.SNP.QUAL," || SOR > ",custom.SNP.SOR," || FS > ",custom.SNP.FS," || MQ < ",custom.SNP.MQ," || MQRankSum < ",custom.SNP.MQRankSum," || ReadPosRankSum < ",custom.SNP.ReadPosRankSum)
    indel.expr=paste0("QD < ",custom.INDEL.QD," || QUAL < ",custom.INDEL.QUAL," || FS > ",custom.INDEL.FS," || ReadPosRankSum < ",custom.INDEL.ReadPosRankSum)
    cmds=c(
      paste(command,"GenotypeGVCFs -R",shQuote(refs[i]),"-V",shQuote(inputs[i]),"--use-new-qual-calculator true -O",shQuote(unfiltered)),
      paste(command,"SelectVariants -R",shQuote(refs[i]),"-V",shQuote(unfiltered),"-O",shQuote(us),"--select-type-to-include SNP"),
      paste(command,"SelectVariants -R",shQuote(refs[i]),"-V",shQuote(unfiltered),"-O",shQuote(ui),"--select-type-to-include INDEL"),
      paste(command,"VariantFiltration -R",shQuote(refs[i]),"-V",shQuote(us),"-O",shQuote(fs),"--filter-name SNP_HARD_FILTER --filter-expression",shQuote(snp.expr)),
      paste(command,"VariantFiltration -R",shQuote(refs[i]),"-V",shQuote(ui),"-O",shQuote(fi),"--filter-name INDEL_HARD_FILTER --filter-expression",shQuote(indel.expr)),
      paste(command,"MergeVcfs -I",shQuote(fs),"-I",shQuote(fi),"-O",shQuote(fg)),
      paste(command,"SelectVariants -R",shQuote(refs[i]),"-V",shQuote(fg),"-O",shQuote(final[1]),"--exclude-filtered true"),
      paste(command,"SelectVariants -R",shQuote(refs[i]),"-V",shQuote(fs),"-O",shQuote(final[2]),"--exclude-filtered true"),
      paste(command,"SelectVariants -R",shQuote(refs[i]),"-V",shQuote(fi),"-O",shQuote(final[3]),"--exclude-filtered true"))
    for(k in seq_along(cmds)) .runCommand(cmds[k],quiet,paste0("genotyping command ",k),stderr.log=log)
    if(!all(file.exists(c(final,indices)))) stop("Final VCF set or indices were not created")
    .markStageComplete(d,"genotypeSamples"); list(success=TRUE)
  },error=function(e) list(success=FALSE,message=conditionMessage(e))),mc.cores=resources$workers)
  .collectWorkers(results,sample.names,"Genotyping"); invisible(sample.names)
}
