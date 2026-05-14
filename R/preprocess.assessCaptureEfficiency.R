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
#' @param bwa.path system path to the directory containing the \code{bwa}
#'   executable; NULL searches the system PATH.
#'
#' @param samtools.path system path to the directory containing the
#'   \code{samtools} executable; NULL searches the system PATH.
#'
#' @param threads number of CPU threads to pass to BWA and samtools.
#'
#' @param mem amount of RAM in GB (currently reserved for future use).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing. Default: \code{FALSE}.
#'
#' @param quiet logical; if TRUE BWA and samtools screen output is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return invisibly; writes per-sample per-target count CSVs to
#'   output.directory and a cross-sample summary to
#'   logs/assessCaptureEfficiency_summary.csv.
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

  # #Debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # input.reads = "processed-reads/cleaned-reads"
  # output.directory = "sample-capture-assessment"
  # target.fasta = "/Users/chutter/Dropbox/Research/0_Github/FrogCap_Files/Probe_Sets/FINAL_marker-seqs_May20-2023.fa"
  # bwa.path = "/Users/chutter/miniconda3/bin"
  # samtools.path = "/Users/chutter/miniconda3/bin"
  # threads = 4
  # mem = 8
  # overwrite = TRUE
  # quiet = TRUE

  # Adds trailing slash to tool paths
  if (is.null(bwa.path) == FALSE) {
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { bwa.path = "" }

  if (is.null(samtools.path) == FALSE) {
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { samtools.path = "" }

  # Quick checks
  if (is.null(input.reads) == TRUE) { stop("Please provide input reads.") }
  if (is.null(target.fasta) == TRUE) { stop("Please provide a target FASTA file.") }
  if (file.exists(target.fasta) == FALSE) { stop("Target FASTA file not found.") }

  # Sets up output directory
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  # Creates log directories
  if (dir.exists("logs/sample_logs") == FALSE) { dir.create("logs/sample_logs", recursive = TRUE) }

  # Builds BWA index of target FASTA once, shared across all samples
  index.path = paste0(output.directory, "/target-index")
  if (dir.exists(index.path) == FALSE) { dir.create(index.path) }
  system(paste0("cp ", target.fasta, " ", index.path, "/targets.fa"))
  system(paste0(bwa.path, "bwa index ", index.path, "/targets.fa"),
         ignore.stdout = quiet, ignore.stderr = quiet)

  # Count total number of target loci in the reference
  n.targets = as.integer(trimws(
    system(paste0("grep -c '^>' ", index.path, "/targets.fa"), intern = TRUE)
  ))

  # Read in sample data
  reads = list.files(input.reads, recursive = TRUE, full.names = TRUE)
  sample.names = list.dirs(input.reads, recursive = FALSE, full.names = FALSE)

  if (length(sample.names) == 0) {
    sample.names = list.files(input.reads, recursive = FALSE, full.names = FALSE)
    sample.names = unique(gsub("_L00.*", "", sample.names))
  }

  # Resumes: skip samples already done
  if (overwrite == FALSE) {
    done.files = list.files(output.directory)
    sample.names = sample.names[!sample.names %in% done.files]
  }

  if (length(sample.names) == 0) { return("No samples remain to analyze.") }

  # Creates the summary log
  summary.data = data.frame(Sample = as.character(),
                            Lane = as.character(),
                            readPairs = as.numeric(),
                            mappedReads = as.numeric(),
                            targetsHit = as.numeric(),
                            totalTargets = as.numeric(),
                            pctTargetsHit = as.numeric(),
                            pctReadsOnTarget = as.numeric(),
                            stringsAsFactors = FALSE)

  for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]

    # Checks the Sample column in case already renamed
    if (length(sample.reads) == 0) { sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }

    sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    # Returns a warning if reads are not found
    if (length(sample.reads) == 0) {
      warning(sample.names[i], " does not have any reads present. Skipping.")
      next
    }#end if

    # Check for empty or near-empty input files (sequencing failures)
    raw.reads = reads[grep(pattern = sample.names[i], x = reads)]
    file.sizes = file.info(raw.reads)$size
    if (any(is.na(file.sizes)) || max(file.sizes, na.rm = TRUE) < 1000) {
      failure.msg = paste0("Sample failed: input read files are empty or near-empty",
                           " (max file size: ", max(file.sizes, na.rm = TRUE), " bytes).",
                           " This indicates a sequencing or library preparation failure.")
      writeLines(failure.msg, paste0("logs/sample_logs/FAILURE_", sample.names[i], ".txt"))
      warning(sample.names[i], " has empty input read files. Skipping.")
      next
    }

    # Creates per-sample output directory
    out.path = paste0(output.directory, "/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path) }

    for (j in 1:length(sample.reads)) {
      #################################################
      ### Part B: map reads to targets with BWA
      #################################################
      lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]

      # Checks in case already renamed
      if (length(lane.reads) == 0) { lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }

      # Returns a warning if reads are not found
      if (length(lane.reads) == 0) {
        warning(sample.reads[j], " does not have any reads present. Skipping.")
        next
      }#end if

      lane.name = gsub(".*/", "", sample.reads[j])

      read1 = lane.reads[grep("_1.f.*|-1.f.*|_R1_.*|-R1_.*|_R1-.*|-R1-.*|READ1.*|_R1.fast.*|-R1.fast.*", lane.reads)]
      read2 = lane.reads[grep("_2.f.*|-2.f.*|_R2_.*|-R2_.*|_R2-.*|-R2-.*|READ2.*|_R2.fast.*|-R2.fast.*", lane.reads)]

      if (length(read1) == 0 || length(read2) == 0) {
        warning(lane.name, " read pairs could not be identified. Skipping.")
        next
      }

      # Maps reads to target sequences
      bam.file = paste0(out.path, "/", lane.name, "_capture.bam")
      system(paste0(bwa.path, "bwa mem -M -t ", threads, " ",
                    index.path, "/targets.fa ",
                    read1[1], " ", read2[1],
                    " | ", samtools.path, "samtools sort -@", threads, " -O BAM",
                    " -o ", bam.file, " -"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      system(paste0(samtools.path, "samtools index ", bam.file),
             ignore.stdout = quiet, ignore.stderr = quiet)

      #################################################
      ### Part C: summarize mapping results
      #################################################
      idx.file = paste0(out.path, "/", lane.name, "_idxstats.txt")
      system(paste0(samtools.path, "samtools idxstats ", bam.file, " > ", idx.file),
             ignore.stdout = quiet, ignore.stderr = quiet)

      idx.data = read.table(idx.file, sep = "\t", header = FALSE,
                            col.names = c("target", "length", "mapped", "unmapped"))
      idx.data = idx.data[idx.data$target != "*", ]

      # Per-target count CSV for detailed inspection
      write.csv(idx.data, file = paste0(out.path, "/", lane.name, "_per-target-counts.csv"),
                row.names = FALSE)

      # Calculates summary statistics
      total.pairs = as.numeric(system(paste0("zcat < ", read1[1], " | echo $((`wc -l`/4))"), intern = TRUE))
      mapped.reads = sum(idx.data$mapped)
      targets.hit = sum(idx.data$mapped > 0)
      pct.targets = round(targets.hit / n.targets * 100, 2)
      pct.on.target = round(mapped.reads / (total.pairs * 2) * 100, 2)

      temp.remove = data.frame(Sample = sample.names[i],
                               Lane = lane.name,
                               readPairs = total.pairs,
                               mappedReads = mapped.reads,
                               targetsHit = targets.hit,
                               totalTargets = n.targets,
                               pctTargetsHit = pct.targets,
                               pctReadsOnTarget = pct.on.target,
                               stringsAsFactors = FALSE)

      summary.data = rbind(summary.data, temp.remove)

      # Removes BAM to save disk space
      system(paste0("rm ", bam.file, " ", bam.file, ".bai"))

      print(paste0(lane.name, " capture assessment complete!"))
    }#end j loop

    print(paste0(sample.names[i], " Completed capture efficiency assessment!"))
  }#end i loop

  write.csv(summary.data, file = "logs/assessCaptureEfficiency_summary.csv", row.names = FALSE)

}#end function
