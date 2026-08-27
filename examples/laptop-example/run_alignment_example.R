#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run with Rscript run_alignment_example.R")
example_directory <- dirname(normalizePath(sub("^--file=", "", script_arg)))

package_source <- Sys.getenv("PHYLOPROCESSR_SOURCE", unset = "")
if (nzchar(package_source)) {
  pkgload::load_all(package_source, quiet = TRUE, export_all = FALSE)
} else {
  suppressPackageStartupMessages(library(PhyloProcessR))
}

source(file.path(example_directory, "configuration.R"), local = TRUE)
setwd(example_directory)
if (!file.exists(file.path(tool_bin, "mafft"))) {
  stop("MAFFT was not found in PHYLOPROCESSR_BIN.")
}
dir.create(file.path(results_directory, "alignments"), recursive = TRUE,
           showWarnings = FALSE)

annotation_base <- file.path(results_directory, "forty-candidate")
annotateTargets(
  assembly.directory = file.path(results_directory, "contigs", "2_reduced-redundancy"),
  target.file = target_markers,
  alignment.contig.name = annotation_base,
  output.directory = file.path(results_directory, "contigs", "4_annotated-contigs"),
  min.match.percent = 60,
  min.match.length = 60,
  min.match.coverage = 50,
  retain.paralogs = FALSE,
  threads = threads,
  memory = memory_gb,
  blast.path = tool_bin,
  cdhit.path = tool_bin,
  overwrite = overwrite,
  quiet = TRUE
)

candidate_directory <- file.path(results_directory, "alignments", "forty-candidate")
alignTargets(
  targets.to.align = paste0(annotation_base, "_to-align.fa"),
  target.file = target_markers,
  output.directory = candidate_directory,
  algorithm = "localpair",
  min.taxa = 4,
  removal.threshold = 0.35,
  threads = threads,
  memory = memory_gb,
  overwrite = overwrite,
  quiet = TRUE,
  mafft.path = tool_bin
)

stats <- read.csv(file.path("logs", "alignTargets_locus_summary.csv"),
                  stringsAsFactors = FALSE)
passed <- stats[stats$Status == "aligned", ]
passed <- passed[order(-passed$N_taxa_aligned, passed$Mean_pct_missing,
                       -passed$Alignment_length, passed$Locus), ]
selected <- head(passed, 20)
selected$Rank <- seq_len(nrow(selected))
selected <- selected[, c("Rank", "Locus", "N_taxa_initial", "N_taxa_aligned",
                         "Alignment_length", "Mean_pct_missing", "Status")]

best_directory <- file.path(results_directory, "alignments", "best-20")
dir.create(best_directory, recursive = TRUE, showWarnings = FALSE)
if (overwrite) unlink(list.files(best_directory, full.names = TRUE))
source_files <- file.path(candidate_directory, paste0(selected$Locus, ".phy"))
if (!all(file.exists(source_files))) stop("A selected alignment file is missing.")
file.copy(source_files, best_directory, overwrite = TRUE)
write.table(selected, file.path(results_directory, "best_20_alignment_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Alignment example completed: ", nrow(passed),
        " candidates passed; ", nrow(selected), " selected.")
