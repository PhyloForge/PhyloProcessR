#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

source("workflow-X2_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
## Workflow X2: Per-sample capture efficiency assessment
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

# Defensive default for configs that pre-date the use.sra setting
if (!exists("use.sra")) use.sra = FALSE

# Set up read source and sample list
if (use.dropbox == TRUE) {

  # Authorize Dropbox connection once at the start
  rdrop2::drop_auth(rdstoken = dropbox.token)

  # Load full sample spreadsheet and get unique sample names.
  # The spreadsheet may have multiple rows per sample (one per lane/file) —
  # dropboxDownload handles multi-lane samples internally, so we loop over
  # unique Sample names, not rows.
  sample.data  = read.csv(sample.file)
  sample.names = unique(sample.data$Sample)
  cat("Found", length(sample.names), "unique samples in", sample.file, "\n")

} else if (use.sra == TRUE) {

  # Read SRA run info and build the same sample names that sraDownload will use,
  # so the loop variable matches what gets downloaded.
  sra.data = read.csv(sra.info.file, stringsAsFactors = FALSE)

  if (!"Run" %in% names(sra.data))
    stop("sra.info.file must contain a 'Run' column with SRR/ERR/DRR accessions.")

  if (!is.null(sra.filter.strategy) && "LibraryStrategy" %in% names(sra.data))
    sra.data = sra.data[sra.data$LibraryStrategy == sra.filter.strategy, ]

  if (!is.null(sra.sample.name.column)) {
    sra.data$sample.name = as.character(sra.data[[sra.sample.name.column]])
  } else if ("ScientificName" %in% names(sra.data)) {
    sci = gsub("[[:space:]]+", "_", trimws(sra.data$ScientificName))
    if ("SampleName" %in% names(sra.data)) {
      sn = trimws(as.character(sra.data$SampleName))
      specimen.id = ifelse(nchar(sn) > 0 & !is.na(sn), sn, sra.data$Run)
      sra.data$sample.name = paste0(sci, "_", specimen.id)
    } else {
      sra.data$sample.name = paste0(sci, "_", sra.data$Run)
    }
  } else {
    sra.data$sample.name = sra.data$Run
  }

  sample.names = unique(sra.data$sample.name)
  cat("Found", length(sample.names), "unique samples in", sra.info.file, "\n")

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

# fastqStats accumulator — written as a rolling CSV after every sample.
# fastpComplete and assessCaptureEfficiency write their own growing CSVs directly.
all.fastq.stats = data.frame()

##################################################################################################
## Main per-sample loop
##################################################################################################

for (i in 1:length(sample.names)) {

  sample.name = sample.names[i]
  cat("\n======================================================\n")
  cat(" Sample", i, "of", length(sample.names), ":", sample.name, "\n")
  cat("======================================================\n")

  assess.dir = paste0("sample-capture-assessment/", sample.name)

  ##############################################################
  ## Resume: if this sample's assessment folder already exists
  ## and contains a per-target CSV, recover its data from the
  ## rolling CSVs and skip re-processing.
  ##############################################################
  if (dir.exists(assess.dir)) {
    done.csvs = list.files(assess.dir, pattern = "_per-target-counts\\.csv$", full.names = TRUE)
    if (length(done.csvs) > 0) {
      cat(" Already processed — reloading saved data and skipping.\n")

      # Reload this sample's fastq stats into the accumulator so the rolling
      # CSV stays complete. fastp and capture stats are read directly from their
      # own growing CSVs at the final merge step — no accumulators needed.
      if (file.exists("logs/X2_fastq-stats_rolling.csv")) {
        tmp = read.csv("logs/X2_fastq-stats_rolling.csv")
        all.fastq.stats = rbind(all.fastq.stats, tmp[tmp$Sample == sample.name, ])
      }
      next
    }
  }

  ##############################################################
  ## Step 1: Obtain reads for this sample
  ##############################################################
  if (use.dropbox == TRUE) {

    # Download from Dropbox — files land flat in raw.dir as:
    #   SampleName_L001_READ1.fastq.gz / SampleName_L001_READ2.fastq.gz
    # Retries up to 2 times if downloaded files fail a gzip integrity check,
    # deleting the corrupt files before each retry so dropboxDownload re-fetches.
    sample.rows = sample.data[sample.data$Sample == sample.name, ]
    temp.csv = tempfile(fileext = ".csv")
    write.csv(sample.rows, temp.csv, row.names = FALSE)

    max.attempts = 3
    download.ok  = FALSE
    last.error   = NULL

    for (attempt in 1:max.attempts) {

      if (attempt > 1) {
        cat(" Retrying download (attempt", attempt, "of", max.attempts, ")...\n")
        Sys.sleep(30)   # longer pause — Dropbox 500 errors usually clear in ~30 s
      }

      # Wrap in tryCatch so Dropbox HTTP errors (500, 503, etc.) don't crash
      # the whole workflow — they're caught here and trigger a retry instead.
      last.error = tryCatch({
        dropboxDownload(sample.spreadsheet = temp.csv,
                        dropbox.directory   = dropbox.directory,
                        dropbox.token       = dropbox.token,
                        output.directory    = raw.dir,
                        overwrite           = TRUE,
                        skip.not.found      = TRUE)
        NULL   # NULL = no error
      }, error = function(e) { conditionMessage(e) })

      if (!is.null(last.error)) {
        cat(" Attempt", attempt, ": download error:", last.error, "\n")
        next   # retry
      }

      # Check files landed
      input.files = list.files(raw.dir, pattern = sample.name, full.names = TRUE)
      input.files = input.files[grep("\\.fastq\\.gz$|\\.fq\\.gz$", input.files)]

      if (length(input.files) == 0) {
        cat(" Attempt", attempt, ": no files found after download.\n")
        last.error = "no files downloaded"
        next
      }

      # gzip integrity check on every downloaded file
      corrupt = sapply(input.files, function(f) {
        system(paste0("gzip -t ", shQuote(f), " 2>/dev/null"), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0
      })

      if (any(corrupt)) {
        cat(" Attempt", attempt, ": corrupted file(s) detected:",
            paste(basename(input.files[corrupt]), collapse = ", "), "\n")
        file.remove(input.files)
        last.error = "corrupted gzip"
      } else {
        download.ok = TRUE
        break
      }

    }#end attempt loop

    unlink(temp.csv)

    if (!download.ok) {
      cat(" WARNING:", sample.name, "skipping after", max.attempts,
          "failed attempts. Last error:", last.error, "\n")
      writeLines(paste0("Download failed after ", max.attempts, " attempts. Last error: ", last.error),
                 paste0("logs/sample_logs/FAILURE_", sample.name, "_download-failed.txt"))
      next
    }

    cat(" Downloaded", length(input.files), "file(s) for", sample.name, "\n")
    input.dir = raw.dir

  } else if (use.sra == TRUE) {

    # Download this sample's single row from ENA via sraDownload.
    # Pass sample.name as a pre-built column so naming is consistent with the
    # sample list built in the setup section above.
    sra.row      = sra.data[sra.data$sample.name == sample.name, ]
    temp.sra.csv = tempfile(fileext = ".csv")
    write.csv(sra.row, temp.sra.csv, row.names = FALSE)

    last.error = tryCatch({
      sraDownload(sra.info.file      = temp.sra.csv,
                  sample.name.column = "sample.name",
                  output.directory   = raw.dir,
                  max.retries        = sra.max.retries,
                  retry.delay        = sra.retry.delay,
                  skip.not.found     = TRUE,
                  overwrite          = FALSE,
                  quiet              = quiet)
      NULL
    }, error = function(e) conditionMessage(e))

    unlink(temp.sra.csv)

    if (!is.null(last.error)) {
      cat(" Download error:", last.error, "\n")
      writeLines(paste0("Download failed. Error: ", last.error),
                 paste0("logs/sample_logs/FAILURE_", sample.name, "_download-failed.txt"))
      next
    }

    input.files = list.files(raw.dir, pattern = sample.name, full.names = TRUE)
    input.files = input.files[grep("\\.fastq\\.gz$|\\.fq\\.gz$", input.files)]

    if (length(input.files) == 0) {
      cat(" WARNING: no files found after SRA download for", sample.name, "— skipping.\n")
      next
    }

    # gzip integrity check
    corrupt = sapply(input.files, function(f) {
      system(paste0("gzip -t ", shQuote(f), " 2>/dev/null"),
             ignore.stdout = TRUE, ignore.stderr = TRUE) != 0
    })
    if (any(corrupt)) {
      cat(" WARNING: corrupted file(s) for", sample.name, "— skipping.\n")
      file.remove(input.files)
      writeLines(paste0("Corrupted gzip: ", paste(basename(input.files[corrupt]), collapse = ", ")),
                 paste0("logs/sample_logs/FAILURE_", sample.name, "_corrupted-file.txt"))
      next
    }

    cat(" Downloaded", length(input.files), "file(s) for", sample.name, "\n")
    input.dir = raw.dir

  } else {

    # Use the local read directory directly — no copying needed.
    # Check whether reads are in a sub-directory or flat.
    sample.subdir = file.path(read.directory, sample.name)
    if (dir.exists(sample.subdir)) {
      input.dir = sample.subdir
    } else {
      input.dir = read.directory
    }

    input.files = list.files(input.dir, pattern = sample.name, full.names = TRUE)
    input.files = input.files[grep("\\.fastq\\.gz$|\\.fq\\.gz$|\\.fastq$|\\.fq$", input.files)]
    if (length(input.files) == 0) {
      cat(" WARNING: no reads found for", sample.name, "in", input.dir, "— skipping.\n")
      next
    }
    cat(" Found", length(input.files), "local file(s) for", sample.name, "\n")

  }

  ##############################################################
  ## Integrity check for local reads: verify gzip files are not
  ## truncated. Dropbox files are already checked during download
  ## with retry above; this catches corruption in local read sets.
  ##############################################################
  if (use.dropbox == FALSE) {
    gz.files = input.files[grep("\\.gz$", input.files)]
    if (length(gz.files) > 0) {
      corrupt = sapply(gz.files, function(f) {
        system(paste0("gzip -t ", shQuote(f), " 2>/dev/null"), ignore.stdout = TRUE, ignore.stderr = TRUE) != 0
      })
      if (any(corrupt)) {
        cat(" WARNING:", sample.name, "has corrupted/truncated local file(s):",
            paste(basename(gz.files[corrupt]), collapse = ", "), "— skipping.\n")
        writeLines(paste0("Corrupted or truncated gzip file(s): ",
                          paste(basename(gz.files[corrupt]), collapse = ", ")),
                   paste0("logs/sample_logs/FAILURE_", sample.name, "_corrupted-file.txt"))
        next
      }
    }
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
  fastpComplete(input.reads       = input.dir,
                output.directory  = cleaned.dir,
                fastp.path        = fastp.path,
                threads           = threads,
                mem               = memory,
                overwrite         = TRUE,
                quiet             = quiet)

  # fastpComplete now appends to its own CSV automatically —
  # no manual accumulation needed here.

  ##############################################################
  ## Step 4a: Barcode identification on cleaned reads
  ## Uses MItoTrawlR::barcodeSampleScan — iterative BBMap read
  ## recruitment + SPAdes/CAP3 assembly stops as soon as the
  ## barcode region is covered, then BLASTs for species ID.
  ##############################################################
  if (run.barcode.scan == TRUE) {
    MItoTrawlR::barcodeSampleScan(
      input.reads    = cleaned.dir,
      output.directory = "barcode-assessment",
      barcode.fasta  = barcode.fasta,
      database.fasta = barcode.database.fasta,
      hits.per.sample = 5,
      per.max.length = barcode.per.max.length,
      min.iterations = barcode.min.iterations,
      max.iterations = barcode.max.iterations,
      min.ref.id     = barcode.min.ref.id,
      bbmap.path     = bbmap.path,
      spades.path    = spades.path,
      cap3.path      = cap3.path,
      blast.path     = blast.path,
      memory         = memory,
      threads        = threads,
      overwrite      = FALSE,
      quiet          = quiet)
  }

  ##############################################################
  ## Step 4b: Assess capture efficiency on cleaned reads
  ## (fastpComplete writes cleaned reads to cleaned.dir/sample/)
  ##############################################################
  assessCaptureEfficiency(input.reads      = cleaned.dir,
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
    raw.files = list.files(raw.dir, pattern = sample.name, full.names = TRUE)
    if (length(raw.files) > 0) {
      file.remove(raw.files)
      cat(" Deleted", length(raw.files), "raw read file(s) for", sample.name, "\n")
    }
  }

  if (use.sra == TRUE && delete.raw.reads == TRUE) {
    raw.files = list.files(raw.dir, pattern = sample.name, full.names = TRUE)
    raw.files = raw.files[grep("\\.fastq\\.gz$|\\.fq\\.gz$", raw.files)]
    if (length(raw.files) > 0) {
      file.remove(raw.files)
      cat(" Deleted", length(raw.files), "raw read file(s) for", sample.name, "\n")
    }
    # Remove the download sentinel so that a re-run (e.g. if assessment failed)
    # triggers a fresh download rather than silently skipping.
    sentinel = file.path(raw.dir, paste0(sample.name, ".fastq.sra_done"))
    if (file.exists(sentinel)) file.remove(sentinel)
  }

  if (delete.cleaned.reads == TRUE) {
    cleaned.sample.dir = paste0(cleaned.dir, "/", sample.name)
    if (dir.exists(cleaned.sample.dir)) {
      system(paste0("rm -rf ", shQuote(cleaned.sample.dir)))
      cat(" Deleted cleaned reads for", sample.name, "\n")
    }
  }

  ##############################################################
  ## Rolling saves — written after every sample so the run can
  ## be safely interrupted and resumed without losing data.
  ##############################################################
  write.csv(all.fastq.stats, "logs/X2_fastq-stats_rolling.csv", row.names = FALSE)
  # fastpComplete and assessCaptureEfficiency write and append their own
  # CSVs directly — logs/fastpComplete_summary.csv and
  # logs/assessCaptureEfficiency_summary.csv grow after every sample

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
# fastpComplete appends to its CSV directly; read it here for the merge.
if (file.exists("logs/fastpComplete_summary.csv")) {
  fp.raw = read.csv("logs/fastpComplete_summary.csv", stringsAsFactors = FALSE)
  fp.agg = aggregate(cbind(startPairs, removePairs, endPairs) ~ Sample,
                     data = fp.raw, FUN = sum)
  fp.agg$pctRemovedByFastp = round(fp.agg$removePairs / fp.agg$startPairs * 100, 2)
  fp.agg = fp.agg[, c("Sample", "startPairs", "removePairs", "endPairs", "pctRemovedByFastp")]
} else {
  fp.agg = data.frame()
}

# --- Capture stats: already aggregated per sample by assessCaptureEfficiency ---
if (file.exists("logs/assessCaptureEfficiency_summary.csv")) {
  cap.agg = read.csv("logs/assessCaptureEfficiency_summary.csv", stringsAsFactors = FALSE)
} else {
  cap.agg = data.frame()
}

# --- Barcode scan: best hit per sample already written by barcodeSampleScan ---
if (run.barcode.scan == TRUE && file.exists("logs/barcodeSampleScan_summary.csv")) {
  bc.agg = read.csv("logs/barcodeSampleScan_summary.csv", stringsAsFactors = FALSE)
} else {
  bc.agg = data.frame()
}

# --- Merge all summaries on Sample ---
merged = fq.agg
if (nrow(fp.agg) > 0)  { merged = merge(merged, fp.agg,  by = "Sample", all = TRUE) }
if (nrow(cap.agg) > 0) { merged = merge(merged, cap.agg, by = "Sample", all = TRUE) }
if (nrow(bc.agg) > 0)  { merged = merge(merged, bc.agg,  by = "Sample", all = TRUE) }

# Reorder columns for readability
keep.cols = intersect(
  c("Sample",
    "Read_Pairs", "MegaBasePairs", "Reads_Per_Million", "Read_Length",
    "startPairs", "removePairs", "endPairs", "pctRemovedByFastp",
    "readPairs", "mappedReads",
    "targetsHit", "totalTargets", "pctTargetsHit", "pctReadsOnTarget",
    "ContigLength", "ContigCount", "BestMatch", "Pident", "AlignLength", "Evalue", "Bitscore"),
  colnames(merged)
)
merged = merged[, keep.cols]

write.csv(merged, "logs/X2_capture-assessment_FINAL.csv", row.names = FALSE)
cat("Final merged summary written to logs/X2_capture-assessment_FINAL.csv\n")
cat("Done! Processed", nrow(merged), "samples.\n")
