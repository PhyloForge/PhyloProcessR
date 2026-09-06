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
#'   logs/assessCaptureEfficiency_summary.csv. Only primary alignments are
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
  # The index is kept between runs. An earlier version copied and indexed the
  # probe set again on every call.
  index.path = paste0(output.directory, "/target-index")
  if (dir.exists(index.path) == FALSE) { dir.create(index.path, recursive = TRUE) }
  target.copy = paste0(index.path, "/targets.fa")

  if (file.exists(target.copy) == FALSE ||
      file.exists(paste0(target.copy, ".bwt")) == FALSE ||
      file.mtime(target.fasta) > file.mtime(target.copy)) {
    file.copy(target.fasta, target.copy, overwrite = TRUE)
    .runCommand(paste0(bwa.command, " index ", shQuote(target.copy)),
                quiet = quiet, task = "bwa index")
  }

  # Count total number of target loci in the reference
  n.targets = as.integer(trimws(
    .runCommandOutput(paste0("grep -c '^>' ", shQuote(target.copy)), task = "target counting")
  ))

  # Read in sample data
  input.reads = sub("/+$", "", input.reads)
  reads = list.files(input.reads, recursive = TRUE, full.names = TRUE)
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

    # Check for empty or near-empty input files (sequencing failures)
    file.sizes = file.info(sample.reads)$size
    file.sizes = file.sizes[is.na(file.sizes) == FALSE]
    if (length(file.sizes) == 0 || max(file.sizes) < 1000) {
      largest.size = if (length(file.sizes) == 0) 0 else max(file.sizes)
      failure.msg = paste0("Sample failed: input read files are empty or near-empty",
                           " (max file size: ", largest.size, " bytes).",
                           " This indicates a sequencing or library preparation failure.")
      writeLines(failure.msg, paste0("logs/sample_logs/FAILURE_", sample.names[i], ".txt"))
      warning(sample.names[i], " has empty input read files. Skipping.")
      next
    }

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

      read1 = lane.reads[grep("_1.f.*|-1.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1.fast.*|-R1.fast.*", basename(lane.reads))]
      read2 = lane.reads[grep("_2.f.*|-2.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2.fast.*|-R2.fast.*", basename(lane.reads))]

      if (length(read1) == 0 || length(read2) == 0) {
        warning(lane.name, " read pairs could not be identified. Skipping.")
        next
      }

      lane.csv = paste0(out.path, "/", lane.name, "_capture-summary.csv")

      # Reuses a finished lane so an interrupted run continues where it stopped
      if (overwrite == FALSE && file.exists(lane.csv) == TRUE) {
        lane.data = rbind(lane.data, read.csv(lane.csv, stringsAsFactors = FALSE))
        print(paste0(lane.name, " is already complete. Skipping."))
        next
      }

      # Maps reads to target sequences. Secondary and supplementary records are
      # dropped here, so a read is counted once. Unmapped records are kept, and
      # they give the read pair total without a second pass over the fastq file.
      bam.file = paste0(out.path, "/", lane.name, "_capture.bam")
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
      idx.file = paste0(out.path, "/", lane.name, "_idxstats.txt")
      .runCommand(paste0(samtools.command, " idxstats ", shQuote(bam.file),
                         " > ", shQuote(idx.file)),
                  quiet = quiet, task = "samtools idxstats", keep.stdout = TRUE)

      idx.data = read.table(idx.file, sep = "\t", header = FALSE,
                            col.names = c("target", "length", "mapped", "unmapped"))
      idx.data = idx.data[idx.data$target != "*", ]

      # Per-target count CSV for detailed inspection
      write.csv(idx.data, file = paste0(out.path, "/", lane.name, "_per-target-counts.csv"),
                row.names = FALSE)

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

      write.csv(temp.remove, file = lane.csv, row.names = FALSE)
      lane.data = rbind(lane.data, temp.remove)

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

  .appendSummary(summary.data, "logs/assessCaptureEfficiency_summary.csv")

  return(invisible(summary.data))
}#end function
