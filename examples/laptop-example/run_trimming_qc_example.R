#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run with Rscript run_trimming_qc_example.R")
example_directory <- dirname(normalizePath(sub("^--file=", "", script_arg)))

package_source <- Sys.getenv("PHYLOPROCESSR_SOURCE", unset = "")
if (nzchar(package_source)) {
  pkgload::load_all(package_source, quiet = TRUE, export_all = FALSE)
} else {
  suppressPackageStartupMessages(library(PhyloProcessR))
}

source(file.path(example_directory, "configuration.R"), local = TRUE)
setwd(example_directory)
if (!file.exists(file.path(tool_bin, "mafft")) ||
    !file.exists(file.path(tool_bin, "trimal"))) {
  stop("MAFFT and TrimAl must both be present in PHYLOPROCESSR_BIN.")
}

input_dir <- file.path(results_directory, "alignments", "best-20")
target_trimmed_dir <- file.path(results_directory, "alignments", "best-20_target-region")
final_dir <- file.path(results_directory, "alignments", "best-20_trimmed")
if (length(list.files(input_dir, pattern = "\\.phy$")) != 20) {
  stop("Run run_alignment_example.R before trimming; 20 input alignments are required.")
}

timings <- data.frame(step = character(), elapsed_seconds = numeric())
time_step <- function(label, expression) {
  elapsed <- system.time(force(expression))[["elapsed"]]
  timings <<- rbind(timings, data.frame(step = label, elapsed_seconds = elapsed))
}

time_step("trim alignments to target regions", trimAlignmentTargets(
  alignment.directory = input_dir,
  alignment.format = "phylip",
  output.directory = target_trimmed_dir,
  target.file = target_markers,
  target.direction = TRUE,
  min.alignment.length = 100,
  min.taxa.alignment = 4,
  threads = threads,
  memory = memory_gb,
  overwrite = overwrite,
  mafft.path = tool_bin
))

time_step("trim and assess target-region alignments", superTrimmer(
  alignment.dir = target_trimmed_dir,
  alignment.format = "phylip",
  output.dir = final_dir,
  TrimAl = TRUE,
  TrimAl.path = tool_bin,
  trim.similarity = TRUE,
  similarity.threshold = 0.4,
  mafft.path = tool_bin,
  trim.external = TRUE,
  min.external.percent = 50,
  trim.coverage = TRUE,
  min.coverage.percent = 35,
  trim.column = TRUE,
  min.column.gap.percent = 50,
  convert.ambiguous.sites = FALSE,
  alignment.assess = TRUE,
  min.coverage.bp = 60,
  min.alignment.length = 100,
  min.taxa.alignment = 4,
  max.alignment.gap.percent = 50,
  threads = threads,
  memory = memory_gb,
  overwrite = overwrite
))
write.csv(timings, file.path(results_directory, "trimming_qc_timings.csv"),
          row.names = FALSE)

best <- read.delim(file.path(results_directory, "best_20_alignment_summary.tsv"),
                   stringsAsFactors = FALSE)
trim <- read.csv(file.path("logs", "best-20_trimmed_trimming_summary.csv"),
                 stringsAsFactors = FALSE)
names(trim)[names(trim) == "Alignment"] <- "Locus"
qc <- merge(best, trim, by = "Locus", sort = FALSE)
qc <- qc[order(qc$Rank), ]
final <- data.frame(
  Rank = qc$Rank,
  Locus = qc$Locus,
  Pass = qc$Pass,
  Untrimmed_taxa = qc$N_taxa_aligned,
  Final_taxa = qc$covSamples,
  Taxa_removed = qc$N_taxa_aligned - qc$covSamples,
  Untrimmed_length_bp = qc$Alignment_length,
  Target_region_length_bp = qc$startLength,
  Final_length_bp = qc$covLength,
  Target_region_retained_pct = round(100 * qc$covLength / qc$startLength, 2),
  Untrimmed_mean_missing_pct = qc$Mean_pct_missing,
  Target_region_gap_pct = round(qc$startPerGaps, 2),
  Final_gap_pct = round(qc$covPerGaps, 2)
)
write.table(final, file.path(results_directory, "final_20_qc_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

files <- list.files(final_dir, pattern = "\\.phy$", full.names = TRUE)
sample_names <- unlist(lapply(files, function(input) names(
  Biostrings::DNAStringSet(Biostrings::readDNAMultipleAlignment(input, format = "phylip"))
)))
sample_summary <- as.data.frame(table(sample_names), stringsAsFactors = FALSE)
names(sample_summary) <- c("Sample", "Final_alignments")
sample_summary <- sample_summary[order(-sample_summary$Final_alignments,
                                       sample_summary$Sample), ]
write.table(sample_summary,
            file.path(results_directory, "final_20_sample_occupancy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Trimming/QC completed: ", sum(final$Pass), " of ", nrow(final),
        " alignments passed; ", sum(final$Taxa_removed), " taxa removed.")
