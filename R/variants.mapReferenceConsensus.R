#' @title mapReferenceConsensus
#'
#' @description Maps per-sample prepared BAM files (from prepareBAM()) against one
#'   shared reference for joint-genotyping workflows, where every sample must use
#'   the same reference, contig names, and coordinates. The reference is built
#'   once by buildReference() into a dataset-owned directory. The reference source
#'   is selected by reference.mode: a majority consensus of the sample alignments
#'   (the default), the capture-target markers FASTA, or a user-supplied FASTA.
#'   Each lane is mapped with the GATK best-practices pipeline
#'   (SamToFastq | bwa mem | MergeBamAlignment | SortSam | MarkDuplicates |
#'   SetNmAndUqTags).
#'
#' @param mapping.directory path to the directory containing per-sample
#'   sub-directories with prepared lane BAM files (all_reads.bam) from
#'   prepareBAM().
#' @param alignment.directory phylip alignment directory used to build the
#'   consensus reference when reference.mode is "consensus".
#' @param samtools.path,bwa.path,gatk4.path tool directories, executable paths,
#'   or NULL to search the system PATH.
#' @param temp.directory path to a GATK JVM temp directory; NULL uses a temporary
#'   directory.
#' @param threads number of CPU threads for BWA and GATK operations.
#' @param memory total JVM heap budget in GB across the two concurrent JVMs in a
#'   mapping pipe.
#' @param overwrite logical; if TRUE the reference and mapped outputs are rebuilt.
#' @param quiet logical; if TRUE tool stdout/stderr is suppressed while logs are
#'   retained.
#' @param reference.path path to the shared reference FASTA to build and map
#'   against. Defaults to the legacy "index/reference.fa"; joint-genotyping
#'   workflows supply a dataset-owned path.
#' @param reference.mode one of "consensus", "target", or "user"; see
#'   buildReference().
#' @param target.file capture-target markers FASTA for reference.mode "target".
#' @param reference.file user-supplied reference FASTA for reference.mode "user".
#' @param sample.names optional retained sample set; NULL discovers every prepared
#'   sample directory.
#'
#' @return invisibly the retained sample names; writes final-mapped-all.bam files
#'   to per-sample lane sub-directories.
#'
#' @export

mapReferenceConsensus = function(mapping.directory = NULL,
                                alignment.directory = NULL,
                                samtools.path = NULL,
                                bwa.path = NULL,
                                gatk4.path = NULL,
                                temp.directory = NULL,
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE,
                                quiet = TRUE,
                                reference.path = "index/reference.fa",
                                reference.mode = c("consensus", "target", "user"),
                                target.file = NULL,
                                reference.file = NULL,
                                sample.names = NULL) {

  reference.mode = match.arg(reference.mode)

  # Quick checks
  if (is.null(mapping.directory) || !dir.exists(mapping.directory)) {
    stop("BAM folder not found.")
  }

  if (is.null(temp.directory)) { temp.directory = tempdir() }
  .ensureDirectory(temp.directory, "temporary directory")
  .ensureDirectory("logs/sample_logs", "sample log directory")

  # Builds and validates the one shared reference before any sample work, so a
  # fully resumed run still confirms the reference is present and unchanged.
  buildReference(reference.path = reference.path,
                 reference.mode = reference.mode,
                 alignment.directory = alignment.directory,
                 target.file = target.file,
                 reference.file = reference.file,
                 samtools.path = samtools.path,
                 bwa.path = bwa.path,
                 gatk4.path = gatk4.path,
                 threads = threads,
                 overwrite = overwrite,
                 quiet = quiet)

  # Discovers prepared samples
  discovered = list.dirs(mapping.directory, recursive = FALSE, full.names = FALSE)
  discovered = discovered[nzchar(discovered)]
  if (is.null(sample.names)) {
    sample.names = discovered
  } else if (any(!sample.names %in% discovered)) {
    stop("Selected sample directories are missing: ",
         paste(setdiff(sample.names, discovered), collapse = ", "))
  }
  if (length(sample.names) == 0) { stop("No prepared samples are available to map.") }

  # Two GATK JVMs run at the same time in a mapping pipe, so the heap budget is
  # split between them.
  resources = .validateResources(threads, memory, 1, simultaneous.jvms = 2)
  gatk = .toolCommand("gatk", gatk4.path)
  bwa = .toolCommand("bwa", bwa.path)
  gatk.command = .gatkCommand(gatk, temp.directory, resources$heap.mb)

  for (sample in sample.names) {
    sample.dir = file.path(mapping.directory, sample)
    lanes = .laneDirectories(sample.dir)
    input.bams = file.path(lanes, "all_reads.bam")
    if (length(lanes) == 0 || any(!file.exists(input.bams))) {
      stop("Incomplete prepared BAM inputs for sample ", sample)
    }

    for (j in seq_along(lanes)) {
      lane.dir = lanes[j]
      output = file.path(lane.dir, "final-mapped-all.bam")
      index = sub("\\.bam$", ".bai", output)
      if (overwrite == FALSE && .stageComplete(lane.dir, "mapReferenceConsensus", c(output, index))) {
        next
      }
      .invalidateStage(lane.dir, "mapReferenceConsensus")
      unlink(c(output, index))

      tmp = file.path(temp.directory, paste0("tmp_", sample, "_", basename(lane.dir)))
      .ensureDirectory(tmp)
      log = file.path("logs/sample_logs", sample, paste0(basename(lane.dir), "_mapping.stderr.log"))
      .ensureDirectory(dirname(log))
      cleaned = file.path(lane.dir, "cleaned_final.bam")
      sorted = file.path(lane.dir, "cleaned_final_sort.bam")
      marked = file.path(lane.dir, "cleaned_final_md.bam")

      # Extracts reads, maps with BWA, and merges the alignment with the unmapped
      # BAM in one streaming pipe.
      pipeline = paste(gatk.command, "SamToFastq -I", shQuote(input.bams[j]),
                       "-FASTQ /dev/stdout -TMP_DIR", shQuote(tmp),
                       "-CLIPPING_ATTRIBUTE XT -CLIPPING_ACTION 2 -INTERLEAVE true -NON_PF true",
                       "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true |",
                       bwa, "mem -M -p -t", threads, shQuote(reference.path), "/dev/stdin |",
                       gatk.command, "MergeBamAlignment -ALIGNED_BAM /dev/stdin -UNMAPPED_BAM",
                       shQuote(input.bams[j]), "-OUTPUT", shQuote(cleaned), "-R", shQuote(reference.path),
                       "-CREATE_INDEX true -ADD_MATE_CIGAR true -CLIP_ADAPTERS false",
                       "-CLIP_OVERLAPPING_READS true -INCLUDE_SECONDARY_ALIGNMENTS true",
                       "-MAX_INSERTIONS_OR_DELETIONS -1 -PRIMARY_ALIGNMENT_STRATEGY MostDistant",
                       "-ATTRIBUTES_TO_RETAIN XS -USE_JDK_DEFLATER true -USE_JDK_INFLATER true",
                       "-TMP_DIR", shQuote(tmp))
      .runPipeline(pipeline, quiet, "read mapping", log)

      # Sort by coordinate for input into MarkDuplicates
      .runCommand(paste(gatk.command, "SortSam -INPUT", shQuote(cleaned), "-OUTPUT", shQuote(sorted),
                        "-CREATE_INDEX true -SORT_ORDER coordinate",
                        "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
                  quiet, "mapping sort", stderr.log = log)
      .runCommand(paste(gatk.command, "MarkDuplicates -INPUT", shQuote(sorted), "-OUTPUT", shQuote(marked),
                        "-CREATE_INDEX true -METRICS_FILE",
                        shQuote(file.path("logs/sample_logs", sample, "duplicate_metrics.txt")),
                        "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
                  quiet, "duplicate marking", stderr.log = log)
      final.pipe = paste(gatk.command, "SortSam -INPUT", shQuote(marked),
                         "-OUTPUT /dev/stdout -SORT_ORDER coordinate |",
                         gatk.command, "SetNmAndUqTags -INPUT /dev/stdin -OUTPUT", shQuote(output),
                         "-CREATE_INDEX true -R", shQuote(reference.path),
                         "-USE_JDK_DEFLATER true -USE_JDK_INFLATER true")
      .runPipeline(final.pipe, quiet, "mapping finalization", log)

      if (!all(file.exists(c(output, index)))) {
        stop("Mapping did not produce BAM and index for ", sample, " ", basename(lane.dir))
      }
      unlink(tmp, recursive = TRUE)
      # GATK -CREATE_INDEX writes <base>.bai, replacing the .bam extension.
      file.remove(c(cleaned, sorted, marked,
                    sub("\\.bam$", ".bai", c(cleaned, sorted, marked))))
      .markStageComplete(lane.dir, "mapReferenceConsensus", paste0("reference=", reference.path))
      if (quiet == FALSE) { message(sample, " ", basename(lane.dir), " completed read mapping to reference.") }
    }#end lane j loop
  }#end sample loop

  invisible(sample.names)
}#end function

# END SCRIPT
