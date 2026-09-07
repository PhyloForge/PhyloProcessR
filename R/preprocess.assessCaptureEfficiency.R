#' @title assessCaptureEfficiency
#'
#' @description Maps cleaned reads back to the target probe/marker sequences
#'   using BWA to estimate sequence capture efficiency for each sample. For
#'   each sample and lane, reports total read pairs, number of reads mapping to
#'   targets, the number of unique target loci with at least one read, and the
#'   percentage of targets recovered and reads on-target. Intended as a quick
#'   QC scan to flag samples with poor enrichment before running the full
#'   assembly pipeline.
#'
#' @param input.reads path to a directory of cleaned reads. Each sample must
#'   occupy its own sub-directory (or be identified by a shared filename
#'   prefix).
#'
#' @param output.directory path to the directory where per-sample mapping files
#'   and per-target count CSVs will be saved. Default:
#'   \code{"sample-capture-assessment"}.
#'
#' @param target.fasta path to the FASTA file of target probe/marker sequences
#'   used for sequence capture.
#'
#' @param bwa.path system path to the directory that contains the \code{bwa}
#'   executable, or the full path to the executable; NULL searches the system
#'   PATH.
#'
#' @param samtools.path system path to the directory that contains the
#'   \code{samtools} executable, or the full path to the executable; NULL
#'   searches the system PATH.
#'
#' @param threads number of CPU threads to pass to BWA and samtools.
#'
#' @param mem amount of RAM in GB, passed to the samtools sort buffer.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing. FALSE resumes and skips only the lanes that
#'   already have a saved result. Default: \code{FALSE}.
#'
#' @param quiet logical; if TRUE BWA and samtools screen output is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return invisibly returns the summary data frame; writes per-sample
#'   per-target count CSVs to output.directory and a cross-sample summary to
#'   logs/sample-capture-assessment_summary.csv. Only primary alignments are
#'   counted, so pctReadsOnTarget cannot go above 100.
#'
#' @export

assessCaptureEfficiency = function(input.reads = NULL,
                                   output.directory = "sample-capture-assessment",
                                   target.fasta = NULL,
                                   bwa.path = NULL,
                                   samtools.path = NULL,
                                   threads = 1,
                                   mem = 8,
                                   overwrite = FALSE,
                                   quiet = TRUE) {

  # Quick checks
  if (is.null(input.reads) == TRUE) { stop("Please provide input reads.") }
  if (file.exists(input.reads) == FALSE) { stop("Input reads not found.") }
  if (is.null(target.fasta) == TRUE) { stop("Please provide a target FASTA file.") }
  if (file.exists(target.fasta) == FALSE) { stop("Target FASTA file not found.") }
  if (length(threads) != 1 || is.numeric(threads) == FALSE ||
      is.finite(threads) == FALSE || threads < 1) {
    stop("threads must be one positive number.")
  }
  if (length(mem) != 1 || is.numeric(mem) == FALSE ||
      is.finite(mem) == FALSE || mem <= 0) {
    stop("mem must be one positive number.")
  }
  if (length(overwrite) != 1 || is.logical(overwrite) == FALSE || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }
  .checkDirectoryOverlap(input.reads, output.directory)

  # Checks that both programs are installed before any sample is processed
  bwa.command = .toolCommand("bwa", bwa.path)
  samtools.command = .toolCommand("samtools", samtools.path)

  # Sets up output directory
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE) { .resetDirectory(output.directory) }
  }#end else

  # Creates log directories
  if (dir.exists("logs/sample_logs") == FALSE) { dir.create("logs/sample_logs", recursive = TRUE) }

  #################################################
  ### Part A: build the target index once
  #################################################
  # The index is used only for this assessment and removed when the function exits.
  index.path = paste0(output.directory, "/target-index")
  on.exit(unlink(index.path, recursive = TRUE), add = TRUE)
  target.copy = paste0(index.path, "/targets.fa")
  target.manifest = paste0(normalizePath(target.fasta), "\t",
                           unname(tools::md5sum(target.fasta)))
  manifest.file = paste0(index.path, "/reference_files.txt")
  index.files = c(target.copy, paste0(target.copy, c(".amb", ".ann", ".bwt", ".pac", ".sa")),
                  manifest.file)
  index.stale = all(file.exists(index.files)) == FALSE
  if (index.stale == FALSE) {
    index.stale = !identical(readLines(manifest.file, warn = FALSE), target.manifest)
  }

  if (index.stale == TRUE) {
    temp.index = tempfile(pattern = "target-index-", tmpdir = output.directory)
    dir.create(temp.index)
    on.exit(unlink(temp.index, recursive = TRUE), add = TRUE)
    temp.target = file.path(temp.index, "targets.fa")
    if (file.copy(target.fasta, temp.target, overwrite = TRUE) == FALSE) {
      stop("Could not copy the target FASTA into the index directory.")
    }
    .runCommand(paste0(bwa.command, " index ", shQuote(temp.target)),
                quiet = quiet, task = "bwa index")
    temp.files = c(temp.target,
                   paste0(temp.target, c(".amb", ".ann", ".bwt", ".pac", ".sa")))
    if (all(file.exists(temp.files)) == FALSE) {
      stop("BWA did not create a complete target index.")
    }
    writeLines(target.manifest, file.path(temp.index, "reference_files.txt"))
    old.index = NULL
    if (dir.exists(index.path) == TRUE) {
      old.index = tempfile(pattern = "target-index-old-", tmpdir = output.directory)
      if (file.rename(index.path, old.index) == FALSE) {
        stop("Could not move the old target index before replacement.")
      }
    }
    if (file.rename(temp.index, index.path) == FALSE) {
      if (is.null(old.index) == FALSE) { file.rename(old.index, index.path) }
      stop("Could not publish the completed target index.")
    }
    if (is.null(old.index) == FALSE) { unlink(old.index, recursive = TRUE) }
  }

  # Count total number of target loci in the reference
  n.targets = as.integer(trimws(
    .runCommandOutput(paste0("grep -c '^>' ", shQuote(target.copy)), task = "target counting")
  ))

  # Read in sample data
  input.reads = sub("/+$", "", input.reads)
  reads = .listFastqFiles(input.reads)
  read.names = .relativePaths(reads, input.reads)
  sample.names = .listSampleNames(input.reads)

  if (length(sample.names) == 0) { return("No samples remain to analyze.") }

  # Creates the per-lane accumulator (collapsed to per-sample before saving)
  lane.data = data.frame(Sample = as.character(),
                         readPairs = as.numeric(),
                         mappedReads = as.numeric(),
                         targetsHit = as.numeric(),
                         totalTargets = as.numeric(),
                         stringsAsFactors = FALSE)

  for (i in seq_along(sample.names)) {
    #################################################
    ### Part B: prepare for loading and checks
    #################################################
    sample.reads = .matchPrefix(reads, read.names, sample.names[i])

    # Returns a warning if reads are not found
    if (length(sample.reads) == 0) {
      warning(sample.names[i], " does not have any reads present. Skipping.")
      next
    }#end if

    # A zero-byte file is incomplete. Small gzip files can contain valid reads.
    failure.file = paste0("logs/sample_logs/FAILURE_", sample.names[i], ".txt")
    file.sizes = file.info(sample.reads)$size
    if (any(is.na(file.sizes)) || any(file.sizes == 0)) {
      writeLines("Sample skipped because at least one input read file is missing or zero bytes.",
                 failure.file)
      warning(sample.names[i], " has a missing or zero-byte input read file. Skipping.")
      next
    }
    if (file.exists(failure.file) == TRUE) { unlink(failure.file) }

    # Creates per-sample output directory
    out.path = paste0(output.directory, "/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path, recursive = TRUE) }

    lane.prefixes = .stripReadSuffix(sample.reads)

    for (j in seq_along(lane.prefixes)) {
      #################################################
      ### Part C: map reads to targets with BWA
      #################################################
      lane.reads = .matchPrefix(reads, reads, lane.prefixes[j])
      lane.name = basename(lane.prefixes[j])

      lane.reads = tryCatch(.orderReadPair(lane.reads), error = function(e) {
        warning(conditionMessage(e))
        return(NULL)
      })
      if (is.null(lane.reads) == TRUE) {
        next
      }
      read1 = lane.reads[1]
      read2 = lane.reads[2]

      lane.csv = paste0(out.path, "/", lane.name, "_capture-summary.csv")
      target.csv = paste0(out.path, "/", lane.name, "_per-target-counts.csv")
      metadata.file = paste0("logs/sample_logs/", sample.names[i], "/",
                             lane.name, "_capture-metadata.csv")
      metadata = .laneMetadata(lane.reads,
                               list(target = target.manifest))

      # Reuses a finished lane so an interrupted run continues where it stopped
      if (overwrite == FALSE && file.exists(lane.csv) == TRUE &&
          file.exists(target.csv) == TRUE &&
          .metadataMatches(metadata.file, metadata) == TRUE) {
        lane.summary = tryCatch(read.csv(lane.csv, stringsAsFactors = FALSE),
                                error = function(e) NULL)
        target.data = tryCatch(read.csv(target.csv, stringsAsFactors = FALSE),
                               error = function(e) NULL)
        if (is.null(lane.summary) == FALSE && is.null(target.data) == FALSE &&
            identical(names(lane.summary), names(lane.data)) &&
            identical(names(target.data), c("target", "length", "mapped", "unmapped"))) {
          lane.data = rbind(lane.data, lane.summary)
          print(paste0(lane.name, " is already complete. Skipping."))
          next
        }
      }
      if (overwrite == FALSE && .metadataConflicts(metadata.file, metadata) == TRUE) {
        stop(lane.name, " was assessed with different inputs or a different target reference. ",
             "Use overwrite = TRUE to replace it.")
      }

      # Maps reads to target sequences. Secondary and supplementary records are
      # dropped here, so a read is counted once. Unmapped records are kept, and
      # they give the read pair total without a second pass over the fastq file.
      bam.file = tempfile(pattern = paste0(lane.name, "-capture-"),
                          tmpdir = out.path, fileext = ".bam")
      on.exit(unlink(c(bam.file, paste0(bam.file, ".bai"))), add = TRUE)
      .runPipeline(paste0(bwa.command, " mem -M -t ", threads, " ",
                          shQuote(target.copy), " ",
                          shQuote(read1[1]), " ", shQuote(read2[1]),
                          " | ", samtools.command, " view -b -F 0x900 - ",
                          " | ", samtools.command, " sort -@ ", threads,
                          " -m ", max(1, floor(mem / max(1, threads))), "G -O BAM",
                          " -o ", shQuote(bam.file), " -"),
                   quiet = quiet, task = "bwa capture mapping")

      .runCommand(paste0(samtools.command, " index ", shQuote(bam.file)),
                  quiet = quiet, task = "samtools index")

      #################################################
      ### Part D: summarize mapping results
      #################################################
      idx.file = tempfile(pattern = paste0(lane.name, "-idxstats-"), tmpdir = out.path)
      on.exit(unlink(idx.file), add = TRUE)
      .runCommand(paste0(samtools.command, " idxstats ", shQuote(bam.file),
                         " > ", shQuote(idx.file)),
                  quiet = quiet, task = "samtools idxstats", keep.stdout = TRUE)

      idx.data = read.table(idx.file, sep = "\t", header = FALSE,
                            col.names = c("target", "length", "mapped", "unmapped"))
      idx.data = idx.data[idx.data$target != "*", ]

      # Per-target count CSV for detailed inspection
      temp.target.csv = tempfile(pattern = paste0(lane.name, "-targets-"),
                                 tmpdir = out.path, fileext = ".csv")
      write.csv(idx.data, file = temp.target.csv, row.names = FALSE)

      # Calculates summary statistics. The first mate of every primary record
      # gives the read pair count.
      total.pairs = as.numeric(.runCommandOutput(
        paste0(samtools.command, " view -c -f 64 ", shQuote(bam.file)),
        task = "samtools read counting"))

      temp.remove = data.frame(Sample = sample.names[i],
                               readPairs = total.pairs,
                               mappedReads = sum(idx.data$mapped),
                               targetsHit = sum(idx.data$mapped > 0),
                               totalTargets = n.targets,
                               stringsAsFactors = FALSE)

      temp.lane.csv = tempfile(pattern = paste0(lane.name, "-summary-"),
                               tmpdir = out.path, fileext = ".csv")
      write.csv(temp.remove, file = temp.lane.csv, row.names = FALSE)
      lane.data = rbind(lane.data, temp.remove)

      .publishFiles(c(temp.target.csv, temp.lane.csv), c(target.csv, lane.csv))
      .writeLaneMetadata(metadata, metadata.file)

      # Removes BAM to save disk space
      unlink(c(bam.file, paste0(bam.file, ".bai")))

      print(paste0(lane.name, " capture assessment complete!"))
    }#end j loop

    print(paste0(sample.names[i], " Completed capture efficiency assessment!"))
  }#end i loop

  # Stops here when every sample was skipped. aggregate() cannot work on an
  # empty table.
  if (nrow(lane.data) == 0) { return("No samples remain to analyze.") }

  # Aggregate lane.data to one row per sample:
  #   readPairs and mappedReads are summed across lanes.
  #   targetsHit uses max across lanes (summing would double-count loci
  #   captured in multiple lanes; per-target CSVs in output.directory can
  #   be used for an exact union if needed).
  sum.agg = aggregate(cbind(readPairs, mappedReads) ~ Sample, data = lane.data, FUN = sum)
  max.agg = aggregate(cbind(targetsHit, totalTargets) ~ Sample, data = lane.data, FUN = max)
  summary.data = merge(sum.agg, max.agg, by = "Sample")
  summary.data$pctTargetsHit    = round(summary.data$targetsHit / summary.data$totalTargets * 100, 2)
  summary.data$pctReadsOnTarget = round(summary.data$mappedReads / (summary.data$readPairs * 2) * 100, 2)

  .appendSummary(summary.data, "logs/sample-capture-assessment_summary.csv")

  return(invisible(summary.data))
}#end function
