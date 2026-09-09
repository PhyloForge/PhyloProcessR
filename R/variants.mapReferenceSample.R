#' Map prepared BAMs to each sample's assembly
#' @param mapping.directory Directory produced by [prepareBAM()].
#' @param assembly.directory Directory of `.fa` or `.fasta` sample assemblies.
#' @param check.assemblies Fail when a prepared sample lacks an assembly.
#' @param samtools.path,bwa.path,gatk4.path Tool directories, executable paths, or NULL.
#' @param threads BWA threads.
#' @param memory Total heap budget in GB (two concurrent JVMs in mapping pipes).
#' @param temp.directory Temporary directory.
#' @param overwrite Rebuild indices and mapped outputs.
#' @param quiet Suppress tool output; diagnostic logs are retained.
#' @return Invisibly returns every retained sample name, including completed ones.
#' @export
mapReferenceSample = function(mapping.directory = NULL, assembly.directory = NULL,
  check.assemblies = TRUE, samtools.path = NULL, bwa.path = NULL, gatk4.path = NULL,
  threads = 1, memory = 1, temp.directory = NULL, overwrite = FALSE, quiet = TRUE) {
  if (is.null(mapping.directory) || !dir.exists(mapping.directory)) stop("BAM folder not found.")
  if (is.null(assembly.directory) || !dir.exists(assembly.directory)) stop("Assembly directory not found.")
  samples = list.dirs(mapping.directory, recursive = FALSE, full.names = FALSE)
  samples = samples[nzchar(samples)]
  assemblies = list.files(assembly.directory, pattern = "\\.(fa|fasta)$", ignore.case = TRUE, full.names = TRUE)
  ids = sub("\\.(fa|fasta)$", "", basename(assemblies), ignore.case = TRUE)
  if (anyDuplicated(ids)) stop("Duplicate assembly sample IDs: ", paste(unique(ids[duplicated(ids)]), collapse = ", "))
  missing = setdiff(samples, ids)
  if (check.assemblies && length(missing)) stop("Missing assemblies for sample(s): ", paste(missing, collapse = ", "))
  retained = intersect(samples, ids)
  if (!length(retained)) stop("No prepared samples have matching assemblies.")
  if (!check.assemblies && length(missing)) message("Skipping samples without assemblies: ", paste(missing, collapse = ", "))
  resources = .validateResources(threads, memory, 1, simultaneous.jvms = 2)
  gatk = .toolCommand("gatk", gatk4.path); bwa = .toolCommand("bwa", bwa.path); samtools = .toolCommand("samtools", samtools.path)
  if (is.null(temp.directory)) temp.directory = tempdir()
  .ensureDirectory(temp.directory, "temporary directory"); .ensureDirectory("logs/sample_logs", "sample log directory")
  gatk.command = .gatkCommand(gatk, temp.directory, resources$heap.mb)

  for (sample in retained) {
    sample.dir = file.path(mapping.directory, sample); index.dir = file.path(sample.dir, "index")
    reference = file.path(index.dir, "reference.fa")
    index.outputs = c(reference, paste0(reference, ".fai"), file.path(index.dir, "reference.dict"), paste0(reference, c(".amb", ".ann", ".bwt", ".pac", ".sa")))
    if (overwrite || !.stageComplete(index.dir, "reference-index", index.outputs)) {
      if (dir.exists(index.dir)) unlink(index.dir, recursive = TRUE)
      .ensureDirectory(index.dir, "reference index directory")
      source = assemblies[match(sample, ids)]
      if (!file.copy(source, reference, overwrite = TRUE)) stop("Could not copy assembly for ", sample)
      log = file.path("logs/sample_logs", sample, "reference-index.stderr.log"); .ensureDirectory(dirname(log))
      .runCommand(paste(bwa, "index -a bwtsw", shQuote(reference)), quiet, "BWA reference indexing", stderr.log = log)
      .runCommand(paste(samtools, "faidx", shQuote(reference)), quiet, "samtools reference indexing", stderr.log = log)
      .runCommand(paste(gatk.command, "CreateSequenceDictionary --REFERENCE", shQuote(reference), "--OUTPUT", shQuote(file.path(index.dir, "reference.dict")), "--USE_JDK_DEFLATER true --USE_JDK_INFLATER true"), quiet, "reference dictionary", stderr.log = log)
      if (!all(file.exists(index.outputs))) stop("Reference indexing did not produce every required file for ", sample)
      .markStageComplete(index.dir, "reference-index", paste0("assembly=", normalizePath(source)))
    }
    lanes = .laneDirectories(sample.dir)
    input.bams = file.path(lanes, "all_reads.bam")
    if (!length(lanes) || any(!file.exists(input.bams))) stop("Incomplete prepared BAM inputs for sample ", sample)
    for (j in seq_along(lanes)) {
      lane.dir = lanes[j]; output = file.path(lane.dir, "final-mapped-all.bam"); index = sub("\\.bam$", ".bai", output)
      if (!overwrite && .stageComplete(lane.dir, "mapReferenceSample", c(output, index))) next
      .invalidateStage(lane.dir, "mapReferenceSample"); file.remove(c(output, index))
      tmp = file.path(temp.directory, paste0("tmp_", sample, "_", basename(lane.dir))); .ensureDirectory(tmp)
      log = file.path("logs/sample_logs", sample, paste0(basename(lane.dir), "_mapping.stderr.log"))
      cleaned = file.path(lane.dir, "cleaned_final.bam"); sorted = file.path(lane.dir, "cleaned_final_sort.bam"); marked = file.path(lane.dir, "cleaned_final_md.bam")
      pipeline = paste(gatk.command, "SamToFastq -I", shQuote(input.bams[j]), "-FASTQ /dev/stdout -TMP_DIR", shQuote(tmp), "-CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true -USE_JDK_DEFLATER true -USE_JDK_INFLATER true |", bwa, "mem -M -p -t", threads, shQuote(reference), "/dev/stdin |", gatk.command, "MergeBamAlignment -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM", shQuote(input.bams[j]), "-OUTPUT", shQuote(cleaned), "-R", shQuote(reference), "-CREATE_INDEX true -ADD_MATE_CIGAR true -CLIP_ADAPTERS false -CLIP_OVERLAPPING_READS true -INCLUDE_SECONDARY_ALIGNMENTS true -MAX_INSERTIONS_OR_DELETIONS -1 -PRIMARY_ALIGNMENT_STRATEGY MostDistant -ATTRIBUTES_TO_RETAIN XS -USE_JDK_DEFLATER true -USE_JDK_INFLATER true -TMP_DIR", shQuote(tmp))
      .runPipeline(pipeline, quiet, "read mapping", log)
      .runCommand(paste(gatk.command, "SortSam -INPUT", shQuote(cleaned), "-OUTPUT", shQuote(sorted), "-CREATE_INDEX true -SORT_ORDER coordinate -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"), quiet, "mapping sort", stderr.log = log)
      .runCommand(paste(gatk.command, "MarkDuplicates -INPUT", shQuote(sorted), "-OUTPUT", shQuote(marked), "-CREATE_INDEX true -METRICS_FILE", shQuote(file.path("logs/sample_logs", sample, "duplicate_metrics.txt")), "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true"), quiet, "duplicate marking", stderr.log = log)
      final.pipe = paste(gatk.command, "SortSam -INPUT", shQuote(marked), "-OUTPUT /dev/stdout -SORT_ORDER coordinate |", gatk.command, "SetNmAndUqTags -INPUT /dev/stdin -OUTPUT", shQuote(output), "-CREATE_INDEX true -R", shQuote(reference), "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true")
      .runPipeline(final.pipe, quiet, "mapping finalization", log)
      if (!all(file.exists(c(output, index)))) stop("Mapping did not produce BAM and index for ", sample, " ", basename(lane.dir))
      unlink(tmp, recursive = TRUE); file.remove(c(cleaned, paste0(cleaned, ".bai"), sorted, paste0(sorted, ".bai"), marked, paste0(marked, ".bai")))
      .markStageComplete(lane.dir, "mapReferenceSample", paste0("reference=", reference))
    }
  }
  invisible(retained)
}
