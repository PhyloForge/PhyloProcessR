#' Prepare paired FASTQ reads as unmapped BAM files
#'
#' @param read.directory Directory containing one directory per sample.
#' @param output.directory Output mapping directory.
#' @param auto.readgroup Derive Illumina read-group fields from the first header.
#' @param samtools.path,bwa.path Retained for compatibility; unused by this step.
#' @param gatk4.path Directory containing `gatk`, its path, or NULL for PATH.
#' @param threads Number of samples processed concurrently.
#' @param memory Total JVM heap budget in GB across concurrent samples.
#' @param temp.directory JVM temporary directory.
#' @param overwrite Recompute owned outputs.
#' @param quiet Suppress tool output while retaining per-lane stderr logs.
#' @return Invisibly returns the sample names.
#' @export
prepareBAM = function(read.directory = NULL, output.directory = "sample-mapping",
                      auto.readgroup = TRUE, samtools.path = NULL, bwa.path = NULL,
                      gatk4.path = NULL, threads = 1, memory = 1,
                      temp.directory = NULL, overwrite = FALSE, quiet = TRUE) {
  if (is.null(read.directory) || !dir.exists(read.directory)) stop("Input read directory not found.")
  samples = .listSampleNames(read.directory)
  if (!length(samples)) stop("No sample directories or FASTQ files were found.")
  resources = .validateResources(threads, memory, length(samples))
  gatk = .toolCommand("gatk", gatk4.path)
  if (is.null(temp.directory)) temp.directory = tempdir()
  .ensureDirectory(temp.directory, "temporary directory")
  if (overwrite && dir.exists(output.directory)) unlink(output.directory, recursive = TRUE)
  .ensureDirectory(output.directory, "output directory")
  .ensureDirectory("logs/sample_logs", "sample log directory")
  gatk.command = .gatkCommand(gatk, temp.directory, resources$heap.mb)

  results = parallel::mclapply(samples, function(sample) {
    tryCatch({
      source.dir = file.path(read.directory, sample)
      if (!dir.exists(source.dir)) source.dir = read.directory
      reads = .listFastqFiles(source.dir, recursive = FALSE)
      if (!length(reads)) stop("no FASTQ files found")
      read.prefixes = vapply(reads, function(x) .stripReadSuffix(x)[1], character(1))
      prefixes = unique(read.prefixes)
      if (identical(source.dir, read.directory))
        prefixes = prefixes[basename(sub("_L[0-9]+$", "", prefixes)) == sample]
      if (!length(prefixes)) stop("no recognizable read pairs found")
      sample.dir = file.path(output.directory, sample)
      .ensureDirectory(sample.dir, "sample output directory")
      for (j in seq_along(prefixes)) {
        lane.reads = reads[read.prefixes == prefixes[j]]
        ordered = tryCatch(.orderReadFiles(lane.reads, allow.third = TRUE), error = function(e)
          stop("lane ", basename(prefixes[j]), ": ", conditionMessage(e)))
        pair = ordered[1:2]
        lane.dir = file.path(sample.dir, paste0("Lane_", j))
        .ensureDirectory(lane.dir, "lane directory")
        output = file.path(lane.dir, "all_reads.bam")
        if (!overwrite && .stageComplete(lane.dir, "prepareBAM", output)) next
        .invalidateStage(lane.dir, "prepareBAM"); file.remove(output)
        log = file.path("logs/sample_logs", sample, paste0("Lane_", j, "_prepareBAM.stderr.log"))
        .ensureDirectory(dirname(log), "sample log directory")
        if (auto.readgroup) {
          con = if (grepl("\\.gz$", pair[1], ignore.case = TRUE)) gzfile(pair[1], "rt") else file(pair[1], "rt")
          header = tryCatch(readLines(con, n = 1, warn = FALSE), finally = close(con))
          fields = if (length(header)) strsplit(sub("^@", "", header), ":", fixed = TRUE)[[1]] else character()
          if (length(header) != 1 || !startsWith(header, "@") || length(fields) < 4 || any(!nzchar(fields[1:4])))
            stop("lane ", basename(prefixes[j]), " has no supported Illumina header; use auto.readgroup = FALSE")
          RGID = paste(fields[3], fields[4], sep = ".")
        } else RGID = paste0("FLOWCELL1.LANE", j)
        RGPU = paste0(RGID, ".", sample)
        fastq.sam = file.path(lane.dir, "fastqsam.bam"); reverted = file.path(lane.dir, "revertsam.bam")
        commands = c(
          paste(gatk.command, "FastqToSam -FASTQ", shQuote(pair[1]), "-FASTQ2", shQuote(pair[2]), "-OUTPUT", shQuote(fastq.sam), "-SAMPLE_NAME", shQuote(sample), "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
          paste(gatk.command, "RevertSam -I", shQuote(fastq.sam), "-O", shQuote(reverted), "-SANITIZE true -MAX_DISCARD_FRACTION 0.005 -ATTRIBUTE_TO_CLEAR XT -ATTRIBUTE_TO_CLEAR XN -ATTRIBUTE_TO_CLEAR AS -ATTRIBUTE_TO_CLEAR OP -SORT_ORDER queryname -RESTORE_ORIGINAL_QUALITIES true -REMOVE_DUPLICATE_INFORMATION true -REMOVE_ALIGNMENT_INFORMATION true -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
          paste(gatk.command, "AddOrReplaceReadGroups -I", shQuote(reverted), "-O", shQuote(output), "-RGSM", shQuote(sample), "-RGPU", shQuote(RGPU), "-RGID", shQuote(RGID), "-RGLB", shQuote(paste0("LIB-", sample)), "-RGPL ILLUMINA -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"))
        for (k in seq_along(commands)) .runCommand(commands[k], quiet, paste0("prepareBAM command ", k), stderr.log = log)
        if (!file.exists(output) || file.info(output)$size == 0) stop("GATK did not create ", output)
        file.remove(c(fastq.sam, reverted))
        .markStageComplete(lane.dir, "prepareBAM", paste0("reads=", paste(pair, collapse = "|")))
      }
      list(success = TRUE)
    }, error = function(e) list(success = FALSE, message = conditionMessage(e)))
  }, mc.cores = resources$workers)
  .collectWorkers(results, samples, "prepareBAM")
  invisible(samples)
}
