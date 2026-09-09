#' Run GATK HaplotypeCaller per sample
#' @param mapping.directory Mapped sample directory.
#' @param output.directory GVCF output directory.
#' @param reference.type `sample` or `consensus`.
#' @param gatk4.path GATK executable directory/path or NULL.
#' @param temp.directory Temporary directory.
#' @param ploidy Positive integer ploidy.
#' @param threads Concurrent samples.
#' @param memory Total JVM heap budget in GB.
#' @param overwrite Recompute outputs.
#' @param quiet Suppress tool output while retaining logs.
#' @param sample.names Optional retained sample set.
#' @return Invisibly returns sample names.
#' @export
haplotypeCaller = function(mapping.directory = NULL, output.directory = "haplotype-caller",
  reference.type = c("sample", "consensus"), gatk4.path = NULL, temp.directory = NULL,
  ploidy = 2, threads = 1, memory = 1, overwrite = FALSE, quiet = TRUE,
  sample.names = NULL) {
  reference.type = match.arg(reference.type)
  if (length(ploidy) != 1 || !is.finite(ploidy) || ploidy < 1 || ploidy != as.integer(ploidy)) stop("ploidy must be a positive integer.")
  if (is.null(mapping.directory) || !dir.exists(mapping.directory)) stop("BAM folder not found.")
  discovered = list.dirs(mapping.directory, recursive = FALSE, full.names = FALSE)
  if (is.null(sample.names)) sample.names = discovered else if (any(!sample.names %in% discovered)) stop("Selected sample directories are missing: ", paste(setdiff(sample.names, discovered), collapse = ", "))
  if (!length(sample.names)) stop("No samples are available to analyze.")
  resources = .validateResources(threads, memory, length(sample.names))
  gatk = .toolCommand("gatk", gatk4.path)
  if (is.null(temp.directory)) temp.directory = tempdir()
  .ensureDirectory(temp.directory); if (overwrite && dir.exists(output.directory)) unlink(output.directory, recursive = TRUE)
  .ensureDirectory(output.directory); .ensureDirectory("logs/sample_logs")
  command = .gatkCommand(gatk, temp.directory, resources$heap.mb)
  results = parallel::mclapply(sample.names, function(sample) tryCatch({
    out.dir = file.path(output.directory, sample); .ensureDirectory(out.dir)
    gvcf = file.path(out.dir, "gatk4-haplotype-caller.g.vcf.gz"); idx = paste0(gvcf, ".tbi")
    if (!overwrite && .stageComplete(out.dir, "haplotypeCaller", c(gvcf, idx))) return(list(success = TRUE))
    .invalidateStage(out.dir, "haplotypeCaller"); file.remove(c(gvcf, idx))
    reference = if (reference.type == "sample") file.path(mapping.directory, sample, "index", "reference.fa") else file.path("index", "reference.fa")
    if (!file.exists(reference)) stop("Reference not found: ", reference)
    lanes = .laneDirectories(file.path(mapping.directory, sample)); bams = file.path(lanes, "final-mapped-all.bam"); bams = bams[file.exists(bams)]
    if (!length(bams)) stop("No completed source lane BAMs found")
    log = file.path("logs/sample_logs", sample, "haplotypeCaller.stderr.log"); .ensureDirectory(dirname(log))
    if (length(bams) > 1) {
      merge.dir = file.path(mapping.directory, sample, "Lane_Merge"); .ensureDirectory(merge.dir)
      merged = file.path(merge.dir, "final-mapped-merge.bam"); sorted = file.path(merge.dir, "final-mapped-sort.bam"); dup = file.path(merge.dir, "final-mapped-dup.bam"); input.bam = file.path(merge.dir, "final-mapped-all.bam")
      inputs = paste(rep("-I", length(bams)), shQuote(bams), collapse = " ")
      .runCommand(paste(command, "MergeSamFiles", inputs, "-O", shQuote(merged), "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true"), quiet, "lane merge", stderr.log = log)
      .runCommand(paste(command, "SortSam -INPUT", shQuote(merged), "-OUTPUT", shQuote(sorted), "-CREATE_INDEX true -SORT_ORDER coordinate -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"), quiet, "merged BAM sort", stderr.log = log)
      .runCommand(paste(command, "MarkDuplicates -INPUT", shQuote(sorted), "-OUTPUT", shQuote(dup), "-CREATE_INDEX true -METRICS_FILE", shQuote(file.path("logs/sample_logs", sample, "duplicate_metrics.txt")), "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true"), quiet, "merged duplicate marking", stderr.log = log)
      .runPipeline(paste(command, "SortSam -INPUT", shQuote(dup), "-OUTPUT /dev/stdout -SORT_ORDER coordinate -USE_JDK_DEFLATER true -USE_JDK_INFLATER true |", command, "SetNmAndUqTags -INPUT /dev/stdin -OUTPUT", shQuote(input.bam), "-CREATE_INDEX true -R", shQuote(reference), "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true"), quiet, "merged BAM finalization", log)
      if (!file.exists(input.bam)) stop("Merged BAM was not created")
      file.remove(c(merged, sorted, dup, paste0(sorted, ".bai"), paste0(dup, ".bai")))
    } else input.bam = bams[1]
    .runCommand(paste(command, "HaplotypeCaller -R", shQuote(reference), "-O", shQuote(gvcf), "-I", shQuote(input.bam), "-ERC GVCF -ploidy", ploidy, "--native-pair-hmm-threads 1 -bamout", shQuote(file.path(out.dir, "gatk4-haplotype-caller.bam"))), quiet, "HaplotypeCaller", stderr.log = log)
    if (!all(file.exists(c(gvcf, idx)))) stop("HaplotypeCaller did not create the GVCF and index")
    .markStageComplete(out.dir, "haplotypeCaller", c(paste0("bam=", input.bam), paste0("ploidy=", ploidy)))
    list(success = TRUE)
  }, error = function(e) list(success = FALSE, message = conditionMessage(e))), mc.cores = resources$workers)
  .collectWorkers(results, sample.names, "HaplotypeCaller")
  invisible(sample.names)
}
