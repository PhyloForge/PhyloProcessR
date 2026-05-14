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
##   1. Download one sample from Dropbox
##   2. Count raw reads (fastqStats)
##   3. Clean reads with fastp (adaptor removal, dedup, length filter)
##   4. Map cleaned reads to probe set and assess capture efficiency
##   5. Delete raw + cleaned reads to free disk space
##   6. Repeat for next sample
##   7. Merge all summaries into a single assessment sheet
##
## Resumable: samples already processed (assessment folder present) are skipped.
## Rolling summary CSVs are saved after each sample so progress is not lost if
## the run is interrupted.
##################################################################################################
##################################################################################################

# Create directory structure
dir.create(processed.reads, showWarnings = FALSE)
dir.create("logs/sample_logs", showWarnings = FALSE, recursive = TRUE)
dir.create("sample-capture-assessment", showWarnings = FALSE)

raw.dir     = paste0(processed.reads, "/raw-reads")
cleaned.dir = paste0(processed.reads, "/cleaned-reads")
dir.create(raw.dir, showWarnings = FALSE)
dir.create(cleaned.dir, showWarnings = FALSE)

# Authorize Dropbox connection
rdrop2::drop_auth(rdstoken = dropbox.token)

# Load full sample spreadsheet
sample.data = read.csv(sample.file)
cat("Found", nrow(sample.data), "samples in", sample.file, "\n")

# Accumulators — built up across iterations and written as rolling CSVs
all.fastq.stats   = data.frame()
all.fastp.stats   = data.frame()
all.capture.stats = data.frame()

##################################################################################################
## Main per-sample loop
##################################################################################################

for (i in 1:nrow(sample.data)) {

  sample.name = sample.data$Sample[i]
  cat("\n======================================================\n")
  cat(" Sample", i, "of", nrow(sample.data), ":", sample.name, "\n")
  cat("======================================================\n")

  assess.dir = paste0("sample-capture-assessment/", sample.name)

  ##############################################################
  ## Resume: if this sample's assessment folder already exists
  ## and contains a per-target CSV, recover its data and skip
  ## the download / process / delete steps.
  ##############################################################
  if (dir.exists(assess.dir)) {
    done.csvs = list.files(assess.dir, pattern = "_per-target-counts\\.csv$", full.names = TRUE)
    if (length(done.csvs) > 0) {
      cat(" Already processed — loading saved summaries and skipping.\n")

      # Reload fastqStats row for this sample if rolling file exists
      if (file.exists("logs/X2_fastq-stats_rolling.csv")) {
        tmp = read.csv("logs/X2_fastq-stats_rolling.csv")
        if (sample.name %in% tmp$Sample) {
          all.fastq.stats = rbind(all.fastq.stats, tmp[tmp$Sample == sample.name, ])
        }
      }
      # Reload fastp row
      if (file.exists("logs/X2_fastp_rolling.csv")) {
        tmp = read.csv("logs/X2_fastp_rolling.csv")
        if (sample.name %in% tmp$Sample) {
          all.fastp.stats = rbind(all.fastp.stats, tmp[tmp$Sample == sample.name, ])
        }
      }
      # Reload capture row
      if (file.exists("logs/X2_capture-assessment_rolling.csv")) {
        tmp = read.csv("logs/X2_capture-assessment_rolling.csv")
        if (sample.name %in% tmp$Sample) {
          all.capture.stats = rbind(all.capture.stats, tmp[tmp$Sample == sample.name, ])
        }
      }
      next
    }
  }

  ##############################################################
  ## Step 1: Download this sample from Dropbox
  ##############################################################
  # Write a temporary single-row spreadsheet so dropboxDownload
  # fetches only this sample instead of the full dataset
  temp.csv = tempfile(fileext = ".csv")
  write.csv(sample.data[i, ], temp.csv, row.names = FALSE)

  dropboxDownload(sample.spreadsheet = temp.csv,
                  dropbox.directory   = dropbox.directory,
                  output.directory    = raw.dir,
                  overwrite           = TRUE,
                  skip.not.found      = TRUE)
  unlink(temp.csv)

  # Verify the download succeeded before continuing
  sample.raw.dir = paste0(raw.dir, "/", sample.name)
  if (!dir.exists(sample.raw.dir) || length(list.files(sample.raw.dir)) == 0) {
    cat(" WARNING:", sample.name, "did not download — skipping.\n")
    next
  }

  ##############################################################
  ## Step 2: FastQ stats on raw reads
  ##############################################################
  fastqStats(read.directory = raw.dir,
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
  ## (adaptor removal, dedup, low-complexity filter, length ≥60)
  ## No decontamination — just enough cleaning to map reliably.
  ##############################################################
  fastpComplete(input.reads       = raw.dir,
                output.directory  = cleaned.dir,
                fastp.path        = fastp.path,
                threads           = threads,
                mem               = memory,
                overwrite         = TRUE,
                quiet             = quiet)

  if (file.exists("logs/fastpComplete_summary.csv")) {
    tmp.fp = read.csv("logs/fastpComplete_summary.csv")
    tmp.fp = tmp.fp[tmp.fp$Sample == sample.name, ]
    all.fastp.stats = rbind(all.fastp.stats, tmp.fp)
  }

  ##############################################################
  ## Step 4: Assess capture efficiency on cleaned reads
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

  if (file.exists("logs/assessCaptureEfficiency_summary.csv")) {
    tmp.cap = read.csv("logs/assessCaptureEfficiency_summary.csv")
    tmp.cap = tmp.cap[tmp.cap$Sample == sample.name, ]
    all.capture.stats = rbind(all.capture.stats, tmp.cap)
  }

  ##############################################################
  ## Step 5: Delete reads for this sample to free disk space
  ##############################################################
  if (delete.raw.reads == TRUE) {
    system(paste0("rm -rf ", shQuote(sample.raw.dir)))
    cat(" Deleted raw reads for", sample.name, "\n")
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
  write.csv(all.fastq.stats,   "logs/X2_fastq-stats_rolling.csv",        row.names = FALSE)
  write.csv(all.fastp.stats,   "logs/X2_fastp_rolling.csv",               row.names = FALSE)
  write.csv(all.capture.stats, "logs/X2_capture-assessment_rolling.csv",  row.names = FALSE)

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
  fq.agg$Read_Length      = all.fastq.stats$Read_Length[1]
  fq.agg$Reads_Per_Million = (fq.agg$Read1_Count + fq.agg$Read2_Count +
                                fq.agg$Read3_Count) / 1000000
} else {
  fq.agg = data.frame()
}

# --- Aggregate fastp stats per sample (sum across lanes) ---
if (nrow(all.fastp.stats) > 0) {
  fp.agg = aggregate(cbind(startPairs, removePairs, endPairs) ~ Sample,
                     data = all.fastp.stats, FUN = sum)
  fp.agg$pctRemovedByFastp = round(fp.agg$removePairs / fp.agg$startPairs * 100, 2)
  # Keep only columns relevant to the merged sheet
  fp.agg = fp.agg[, c("Sample", "startPairs", "removePairs", "endPairs", "pctRemovedByFastp")]
} else {
  fp.agg = data.frame()
}

# --- Aggregate capture stats per sample ---
# readPairs and mappedReads are summed across lanes.
# targetsHit is taken as the MAX across lanes — summing would double-count loci
# captured in multiple lanes. The per-target CSVs in sample-capture-assessment/
# can be inspected to compute exact union-of-lanes coverage if needed.
if (nrow(all.capture.stats) > 0) {
  cap.sum = aggregate(cbind(readPairs, mappedReads) ~ Sample,
                      data = all.capture.stats, FUN = sum)
  cap.max = aggregate(cbind(targetsHit, totalTargets) ~ Sample,
                      data = all.capture.stats, FUN = max)
  cap.agg = merge(cap.sum, cap.max, by = "Sample")
  cap.agg$pctTargetsHit    = round(cap.agg$targetsHit / cap.agg$totalTargets * 100, 2)
  cap.agg$pctReadsOnTarget = round(cap.agg$mappedReads / (cap.agg$readPairs * 2) * 100, 2)
  cap.agg = cap.agg[, c("Sample", "readPairs", "mappedReads",
                         "targetsHit", "totalTargets",
                         "pctTargetsHit", "pctReadsOnTarget")]
} else {
  cap.agg = data.frame()
}

# --- Merge the three summaries on Sample ---
merged = fq.agg
if (nrow(fp.agg) > 0)  { merged = merge(merged, fp.agg,  by = "Sample", all = TRUE) }
if (nrow(cap.agg) > 0) { merged = merge(merged, cap.agg, by = "Sample", all = TRUE) }

# Reorder columns for readability
keep.cols = intersect(
  c("Sample",
    "Read_Pairs", "MegaBasePairs", "Reads_Per_Million", "Read_Length",
    "startPairs", "removePairs", "endPairs", "pctRemovedByFastp",
    "readPairs", "mappedReads",
    "targetsHit", "totalTargets", "pctTargetsHit", "pctReadsOnTarget"),
  colnames(merged)
)
merged = merged[, keep.cols]

write.csv(merged, "logs/X2_capture-assessment_FINAL.csv", row.names = FALSE)
cat("Final merged summary written to logs/X2_capture-assessment_FINAL.csv\n")
cat("Done! Processed", nrow(merged), "samples.\n")
