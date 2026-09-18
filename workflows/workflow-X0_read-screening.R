source("workflow-X0_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)
if (!file.exists(target.fasta)) {
  stop("Target FASTA file not found: ", target.fasta)
}
target.md5 = unname(tools::md5sum(target.fasta))

##################################################################################################
##################################################################################################
## Workflow X0: Rapid per-sample read screening
##
## Processes one sample at a time to avoid storing the full dataset on disk:
##   1. Download one sample from Dropbox (files land flat in raw-reads/ as
##      SampleName_L001_READ1.fastq.gz — no per-sample subdirectory)
##   2. Count raw reads (fastqStats)
##   3. Clean reads with fastp (adaptor removal, dedup, length filter)
##   4. Map cleaned reads to probe set and assess capture efficiency
##   5. Delete raw + cleaned reads to free disk space
##   6. Repeat for next sample
##   7. Merge all summaries into a single assessment sheet
##
## Resumable: samples already processed (assessment folder present) are skipped.
## Rolling summary CSVs are written after every sample so progress is not lost
## if the run is interrupted.
##################################################################################################
##################################################################################################

# Create directory structure
dir.create(processed.reads, showWarnings = FALSE)
dir.create("logs/sample_logs", showWarnings = FALSE, recursive = TRUE)
dir.create("sample-capture-assessment", showWarnings = FALSE)
if (run.barcode.scan == TRUE) { dir.create("barcode-assessment", showWarnings = FALSE) }

raw.dir     = paste0(processed.reads, "/raw-reads")
cleaned.dir = paste0(processed.reads, "/cleaned-reads")
dir.create(raw.dir, showWarnings = FALSE)
dir.create(cleaned.dir, showWarnings = FALSE)

# Set up read source and sample list
if (use.dropbox == TRUE) {
  # Load full sample spreadsheet and get unique sample names.
  # The spreadsheet may have multiple rows per sample (one per lane/file) —
  # dropboxDownload handles multi-lane samples internally, so we loop over
  # unique Sample names, not rows.
  sample.data  = read.csv(sample.file, stringsAsFactors = FALSE)
  if (!all(c("File", "Sample") %in% names(sample.data))) {
    stop("sample.file must contain File and Sample columns: ", sample.file)
  }
  sample.names = unique(sample.data$Sample)
  cat("Found", length(sample.names), "unique samples in", sample.file, "\n")

} else {

  # Discover samples from the local read directory.
  # Supports both sub-directory-per-sample and flat file layouts.
  sample.names = list.dirs(read.directory, recursive = FALSE, full.names = FALSE)
  if (length(sample.names) == 0) {
    local.files  = list.files(read.directory, recursive = FALSE, full.names = FALSE)
    sample.names = unique(gsub("_L00.*|_R[12].*|_READ[12].*", "", local.files))
    sample.names = sample.names[grep("\\.fastq|\\.fq", sample.names, invert = TRUE)]
  }
  cat("Found", length(sample.names), "unique samples in", read.directory, "\n")

}

sample.names = sample.names[!is.na(sample.names) & nzchar(sample.names)]
if (length(sample.names) == 0) {
  stop("No sample names were found in the configured read source.")
}

# Record the requested samples even when a sample fails before a summary exists.
sample.status = data.frame(Sample = sample.names,
                           Status = rep("pending", length(sample.names)),
                           stringsAsFactors = FALSE)

# fastqStats accumulator — written as a rolling CSV after every sample.
# fastpClean and assessCaptureEfficiency write their own growing CSVs directly.
all.fastq.stats = data.frame()

##################################################################################################
## Main per-sample loop
##################################################################################################

for (i in seq_along(sample.names)) {

  sample.name = sample.names[i]
  cat("\n======================================================\n")
  cat(" Sample", i, "of", length(sample.names), ":", sample.name, "\n")
  cat("======================================================\n")

  assess.dir = paste0("sample-capture-assessment/", sample.name)
  completion.file = file.path("logs", "sample_logs", sample.name,
                              paste0(sample.name, "_X0-complete.csv"))

  ##############################################################
  ## Resume: if this sample's assessment folder already exists
  ## and contains a per-target CSV, recover its data from the
  ## rolling CSVs and skip re-processing.
  ##############################################################
  completion.matches = FALSE
  if (file.exists(completion.file)) {
    completion = tryCatch(read.csv(completion.file, stringsAsFactors = FALSE),
                          error = function(e) NULL)
    completion.matches = !is.null(completion) &&
      identical(as.character(completion$TargetMD5), target.md5)
  }

  if (dir.exists(assess.dir) && completion.matches) {
    done.csvs = list.files(assess.dir, pattern = "_per-target-counts\\.csv$", full.names = TRUE)
    if (length(done.csvs) > 0) {
      cat(" Already complete — reloading saved data and skipping.\n")
      sample.status$Status[sample.status$Sample == sample.name] = "complete"

      # Reload this sample's fastq stats into the accumulator so the rolling
      # CSV stays complete. fastp and capture stats are read directly from their
      # own growing CSVs at the final merge step — no accumulators needed.
      if (file.exists("logs/X0_fastq-stats_rolling.csv")) {
        tmp = read.csv("logs/X0_fastq-stats_rolling.csv")
        all.fastq.stats = rbind(all.fastq.stats, tmp[tmp$Sample == sample.name, ])
      }
      next
    }
  }

  # A changed target reference invalidates this sample's capture assessment.
  if (dir.exists(assess.dir) && !completion.matches) {
    unlink(assess.dir, recursive = TRUE)
    sample.log.dir = file.path("logs", "sample_logs", sample.name)
    capture.metadata = list.files(sample.log.dir,
                                  pattern = "_capture-metadata\\.csv$",
                                  full.names = TRUE)
    if (length(capture.metadata) > 0) unlink(capture.metadata)
  }

  ##############################################################
  ## Step 1: Obtain reads for this sample
  ##############################################################
  if (use.dropbox == TRUE) {

    # Keep each download in its own directory. A stopped run can leave reads
    # behind, and those reads must not enter the next sample's fastp run.
    sample.raw.dir = file.path(raw.dir, sample.name)
    dir.create(sample.raw.dir, showWarnings = FALSE, recursive = TRUE)

    sample.rows = sample.data[sample.data$Sample == sample.name, ]
    temp.csv = tempfile(fileext = ".csv")
    write.csv(sample.rows, temp.csv, row.names = FALSE)

    dropboxDownload(sample.spreadsheet = temp.csv,
                    dropbox.directory   = dropbox.directory,
                    dropbox.token       = dropbox.token,
                    output.directory    = sample.raw.dir,
                    overwrite           = FALSE,
                    skip.not.found      = TRUE)
    unlink(temp.csv)

    # Verify download succeeded
    input.files = list.files(sample.raw.dir, pattern = sample.name, full.names = TRUE)
    input.files = input.files[grep("\\.fastq\\.gz$|\\.fq\\.gz$", input.files)]
    if (length(input.files) == 0) {
      cat(" WARNING:", sample.name, "did not download — skipping.\n")
      sample.status$Status[sample.status$Sample == sample.name] = "failed: download"
      next
    }
    cat(" Downloaded", length(input.files), "file(s) for", sample.name, "\n")
    input.dir = sample.raw.dir

  } else {

    # Use an existing sample directory when present. For a flat input layout,
    # make a directory of links so package functions see only this sample.
    sample.subdir = file.path(read.directory, sample.name)
    local.staging = FALSE
    if (dir.exists(sample.subdir)) {
      input.dir = sample.subdir
      input.files = list.files(input.dir, full.names = TRUE)
    } else {
      local.files = list.files(read.directory, recursive = FALSE, full.names = TRUE)
      local.names = basename(local.files)
      sample.prefixes = c(paste0(sample.name, "_"), paste0(sample.name, "-"),
                          paste0(sample.name, "."))
      is.sample = Reduce(`|`, lapply(sample.prefixes, function(prefix) {
        startsWith(local.names, prefix)
      }))
      input.files = local.files[is.sample]
      input.files = input.files[grep("\\.fastq(\\.gz)?$|\\.fq(\\.gz)?$",
                                     input.files, ignore.case = TRUE)]
      input.dir = file.path(raw.dir, sample.name)
      if (dir.exists(input.dir)) unlink(input.dir, recursive = TRUE)
      dir.create(input.dir, recursive = TRUE, showWarnings = FALSE)
      linked = file.symlink(input.files, file.path(input.dir, basename(input.files)))
      if (length(linked) > 0 && !all(linked)) {
        stop("Could not stage the local reads for ", sample.name, ".")
      }
      input.files = list.files(input.dir, full.names = TRUE)
      local.staging = TRUE
    }

    input.files = input.files[grep("\\.fastq\\.gz$|\\.fq\\.gz$|\\.fastq$|\\.fq$", input.files)]
    if (length(input.files) == 0) {
      cat(" WARNING: no reads found for", sample.name, "in", input.dir, "— skipping.\n")
      sample.status$Status[sample.status$Sample == sample.name] = "failed: no reads"
      next
    }
    cat(" Found", length(input.files), "local file(s) for", sample.name, "\n")

  }

  ##############################################################
  ## Integrity check: verify gzip files are not truncated.
  ## A corrupted download produces "unexpected end of file" errors
  ## that cause fastp to hang or crash mid-run.
  ##############################################################
  gz.files = input.files[grep("\\.gz$", input.files)]
  if (length(gz.files) > 0) {
    corrupt = sapply(gz.files, function(f) {
      system(paste0("gzip -t ", shQuote(f), " 2>/dev/null"), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0
    })
    if (any(corrupt)) {
      cat(" WARNING:", sample.name, "has corrupted/truncated file(s):",
          paste(basename(gz.files[corrupt]), collapse = ", "), "— skipping.\n")
      writeLines(paste0("Corrupted or truncated gzip file(s): ",
                        paste(basename(gz.files[corrupt]), collapse = ", ")),
                 paste0("logs/sample_logs/FAILURE_", sample.name, "_corrupted-download.txt"))
      # Remove the bad files so a future re-run re-downloads them
      if (use.dropbox == TRUE) { file.remove(input.files) }
      sample.status$Status[sample.status$Sample == sample.name] = "failed: corrupt reads"
      next
    }
  }

  # A sample is one screening unit. Do not estimate success from only the
  # complete lanes when another lane has a missing mate.
  strip.read.suffix = getFromNamespace(".stripReadSuffix", "PhyloProcessR")
  match.prefix = getFromNamespace(".matchPrefix", "PhyloProcessR")
  order.read.pair = getFromNamespace(".orderReadPair", "PhyloProcessR")
  lane.prefixes = strip.read.suffix(input.files)
  complete.lanes = vapply(lane.prefixes, function(lane.prefix) {
    lane.files = match.prefix(input.files, input.files, lane.prefix)
    tryCatch(length(order.read.pair(lane.files)) == 2,
             error = function(e) FALSE)
  }, logical(1))
  if (length(complete.lanes) == 0 || !all(complete.lanes)) {
    cat(" WARNING:", sample.name, "has an incomplete read pair — skipping.\n")
    sample.status$Status[sample.status$Sample == sample.name] =
      "failed: incomplete read pair"
    if (exists("local.staging") && isTRUE(local.staging) && dir.exists(input.dir)) {
      unlink(input.dir, recursive = TRUE)
    }
    next
  }

  ##############################################################
  ## Step 2: FastQ stats on raw reads
  ##############################################################
  fastqStats(read.directory = input.dir,
             output.name    = "fastq-stats-temp",
             read.length    = read.length,
             threads        = threads,
             mem            = memory,
             overwrite      = TRUE)

  if (file.exists("fastq-stats-temp.csv")) {
    tmp.fq = read.csv("fastq-stats-temp.csv")
    tmp.fq = tmp.fq[tmp.fq$Sample == sample.name, ]
    all.fastq.stats = rbind(all.fastq.stats, tmp.fq)
    unlink("fastq-stats-temp.csv")
  }

  ##############################################################
  ## Step 3: Clean reads with fastp
  ## (adaptor removal, dedup, low-complexity filter, length >=60)
  ## No decontamination — just enough cleaning to map reliably.
  ##############################################################
  fastpClean(input.reads           = input.dir,
             output.directory      = cleaned.dir,
             remove.adaptors       = fastp.remove.adaptors,
             remove.duplicate.reads = fastp.remove.duplicate.reads,
             error.correction      = fastp.error.correction,
             quality.trim.reads    = fastp.quality.trim.reads,
             quality.filter        = fastp.quality.filter,
             low.complexity.filter = fastp.low.complexity.filter,
             trim.poly.x           = fastp.trim.poly.x,
             min.read.length       = fastp.min.read.length,
             fastp.path            = fastp.path,
             threads               = threads,
             mem                   = memory,
             overwrite             = FALSE,
             quiet                 = quiet)

  # fastpClean appends its summary to logs/fastp_summary.csv.

  cleaned.sample.dir = file.path(cleaned.dir, sample.name)
  cleaned.files = list.files(cleaned.sample.dir,
                             pattern = "\\.fastq(\\.gz)?$|\\.fq(\\.gz)?$",
                             full.names = TRUE, ignore.case = TRUE)
  if (length(cleaned.files) < 2) {
    cat(" WARNING:", sample.name, "did not produce a cleaned read pair — skipping.\n")
    sample.status$Status[sample.status$Sample == sample.name] = "failed: incomplete read pair"
    next
  }

  ##############################################################
  ## Step 4a: Barcode identification on cleaned reads
  ## Maps reads to a barcode reference (e.g. 16S, COI), assembles
  ## on-target reads with SPAdes, and identifies the best species
  ## match via BLAST. Runs on the same cleaned reads as the capture
  ## assessment so no extra cleaning step is needed.
  ##############################################################
  if (run.barcode.scan == TRUE) {
    if (!requireNamespace("MItoTrawlR", quietly = TRUE)) {
      stop("run.barcode.scan = TRUE requires the MItoTrawlR package.")
    }
    MItoTrawlR::barcodeSampleScan(
      input.reads = cleaned.sample.dir,
      output.directory = "barcode-assessment",
      barcode.fasta = barcode.fasta,
      database.fasta = barcode.database.fasta,
      hits.per.sample = barcode.hits.per.sample,
      per.max.length = barcode.per.max.length,
      min.iterations = barcode.min.iterations,
      max.iterations = barcode.max.iterations,
      min.ref.id = barcode.min.ref.id,
      bbmap.path = bbmap.path,
      spades.path = spades.path,
      cap3.path = cap3.path,
      blast.path = blast.path,
      memory = memory,
      threads = threads,
      overwrite = FALSE,
      quiet = quiet
    )
  }

  ##############################################################
  ## Step 4b: Assess capture efficiency on cleaned reads
  ## Assess only the current sample's cleaned-read directory.
  ##############################################################
  assessCaptureEfficiency(input.reads      = cleaned.sample.dir,
                          output.directory = "sample-capture-assessment",
                          target.fasta     = target.fasta,
                          bwa.path         = bwa.path,
                          samtools.path    = samtools.path,
                          threads          = threads,
                          mem              = memory,
                          overwrite        = FALSE,   # accumulate per-target CSVs across samples
                          quiet            = quiet)

  # assessCaptureEfficiency now appends to its own CSV automatically —
  # no manual accumulation needed here.

  ##############################################################
  ## Step 5: Delete reads for this sample to free disk space
  ## Raw reads are only deleted when using Dropbox (local reads
  ## are never touched regardless of delete.raw.reads).
  ##############################################################
  if (use.dropbox == TRUE && delete.raw.reads == TRUE) {
    raw.files = list.files(input.dir, full.names = TRUE)
    if (dir.exists(input.dir)) {
      unlink(input.dir, recursive = TRUE)
      cat(" Deleted", length(raw.files), "raw read file(s) for", sample.name, "\n")
    }
  }

  dir.create(dirname(completion.file), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(Sample = sample.name,
                       TargetMD5 = target.md5,
                       stringsAsFactors = FALSE),
            completion.file, row.names = FALSE)
  sample.status$Status[sample.status$Sample == sample.name] = "complete"

  if (delete.cleaned.reads == TRUE) {
    if (dir.exists(cleaned.sample.dir)) {
      unlink(cleaned.sample.dir, recursive = TRUE)
      cat(" Deleted cleaned reads for", sample.name, "\n")
    }
  }

  if (exists("local.staging") && isTRUE(local.staging) && dir.exists(input.dir)) {
    unlink(input.dir, recursive = TRUE)
  }

  ##############################################################
  ## Rolling saves — written after every sample so the run can
  ## be safely interrupted and resumed without losing data.
  ##############################################################
  write.csv(all.fastq.stats, "logs/X0_fastq-stats_rolling.csv", row.names = FALSE)
  # fastpClean and assessCaptureEfficiency append their own summary CSVs.

  cat(" Sample", sample.name, "complete!\n")

}#end sample loop

##################################################################################################
## Step 6: Merge all summaries into one assessment sheet
##################################################################################################
cat("\nMerging summaries...\n")

# --- Aggregate fastqStats per sample (sum across lanes) ---
if (nrow(all.fastq.stats) > 0) {
  fq.agg = aggregate(cbind(Read1_Count, Read2_Count, Read3_Count,
                            Total_Reads, Read_Pairs, MegaBasePairs) ~ Sample,
                     data = all.fastq.stats, FUN = sum)
  fq.agg$Read_Length       = all.fastq.stats$Read_Length[1]
  fq.agg$Reads_Per_Million = (fq.agg$Read1_Count + fq.agg$Read2_Count +
                                fq.agg$Read3_Count) / 1000000
} else {
  fq.agg = data.frame()
}

# --- Aggregate fastp stats per sample (sum across lanes) ---
# fastpClean appends to its CSV directly; read it here for the merge.
if (file.exists("logs/fastp_summary.csv")) {
  fp.raw = read.csv("logs/fastp_summary.csv", stringsAsFactors = FALSE)
  fp.agg = aggregate(cbind(startPairs, removePairs, endPairs) ~ Sample,
                     data = fp.raw, FUN = sum)
  fp.agg$pctRemovedByFastp = round(fp.agg$removePairs / fp.agg$startPairs * 100, 2)
  fp.agg = fp.agg[, c("Sample", "startPairs", "removePairs", "endPairs", "pctRemovedByFastp")]
} else {
  fp.agg = data.frame()
}

# --- Capture stats: already aggregated per sample by assessCaptureEfficiency ---
if (file.exists("logs/sample-capture-assessment_summary.csv")) {
  cap.agg = read.csv("logs/sample-capture-assessment_summary.csv", stringsAsFactors = FALSE)
} else {
  cap.agg = data.frame()
}

# --- Barcode scan: best hit per sample already written by barcodeSampleScan ---
if (run.barcode.scan == TRUE && file.exists("logs/barcodeSampleScan_summary.csv")) {
  bc.agg = read.csv("logs/barcodeSampleScan_summary.csv", stringsAsFactors = FALSE)
} else {
  bc.agg = data.frame()
}

merged.tables = list(sample.status, fq.agg, fp.agg, cap.agg, bc.agg)
merged.tables = merged.tables[vapply(merged.tables, nrow, integer(1)) > 0]
merged = Reduce(function(x, y) merge(x, y, by = "Sample", all = TRUE), merged.tables)

# Reorder columns for readability.
keep.cols = intersect(
  c("Sample", "Status",
    "Read_Pairs", "MegaBasePairs", "Reads_Per_Million", "Read_Length",
    "startPairs", "removePairs", "endPairs", "pctRemovedByFastp",
    "readPairs", "mappedReads",
    "targetsHit", "totalTargets", "pctTargetsHit", "pctReadsOnTarget",
    "MappedReads", "ContigLength", "BestMatch", "Pident", "AlignLength", "Evalue", "Bitscore"),
  colnames(merged)
)
merged = merged[, keep.cols, drop = FALSE]

write.csv(merged, "logs/X0_read-screening_FINAL.csv", row.names = FALSE)
cat("Final merged summary written to logs/X0_read-screening_FINAL.csv\n")
cat("Done! Processed", nrow(merged), "samples.\n")
