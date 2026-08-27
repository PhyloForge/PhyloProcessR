#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Run this example with Rscript run_laptop_example.R")
}

example_directory <- dirname(
  normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
)

package_source <- Sys.getenv("PHYLOPROCESSR_SOURCE", unset = "")
if (nzchar(package_source)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required when PHYLOPROCESSR_SOURCE is used.")
  }
  pkgload::load_all(package_source, quiet = TRUE, export_all = FALSE)
} else {
  suppressPackageStartupMessages(library(PhyloProcessR))
}

source(file.path(example_directory, "configuration.R"), local = TRUE)
setwd(example_directory)

dir.create(results_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(results_directory, "processed-reads"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(file.path(results_directory, "contigs"), recursive = TRUE,
           showWarnings = FALSE)

timings <- data.frame(
  step = character(),
  elapsed_seconds = numeric(),
  stringsAsFactors = FALSE
)

run_step <- function(step_name, expression) {
  message("\n--- ", step_name, " ---")
  elapsed <- system.time(force(expression))[["elapsed"]]
  timings <<- rbind(
    timings,
    data.frame(step = step_name, elapsed_seconds = elapsed)
  )
  invisible(elapsed)
}

run_step(
  "summarize raw reads",
  fastqStats(
    read.directory = raw_reads,
    output.name = file.path(results_directory, "fastq-stats"),
    read.length = 101,
    threads = threads,
    mem = memory_gb,
    overwrite = overwrite
  )
)

cleaned_reads <- file.path(results_directory, "processed-reads", "cleaned-reads")
run_step(
  "clean paired reads with fastp",
  fastpComplete(
    input.reads = raw_reads,
    output.directory = cleaned_reads,
    fastp.path = tool_bin,
    threads = threads,
    mem = memory_gb,
    overwrite = overwrite,
    quiet = TRUE
  )
)

draft_contigs <- file.path(results_directory, "contigs", "1_draft-contigs")
run_step(
  "assemble cleaned reads with SPAdes",
  assembleSpades(
    input.reads = cleaned_reads,
    output.directory = file.path(results_directory, "spades-work"),
    assembly.directory = draft_contigs,
    spades.path = tool_bin,
    mismatch.corrector = spades_mismatch_corrector,
    isolate = FALSE,
    kmer.values = spades_kmers,
    threads = threads,
    memory = memory_gb,
    overwrite = overwrite,
    save.corrected.reads = FALSE,
    clean.up.spades = TRUE,
    quiet = TRUE
  )
)

reduced_contigs <- file.path(results_directory, "contigs", "2_reduced-redundancy")
run_step(
  "reduce redundant contigs",
  reduceRedundancy(
    assembly.directory = draft_contigs,
    output.directory = reduced_contigs,
    similarity = reduce_redundancy_similarity,
    cdhit.path = tool_bin,
    memory = memory_gb,
    threads = threads,
    overwrite = overwrite,
    quiet = TRUE
  )
)

target_contigs <- file.path(results_directory, "contigs", "3_target-contigs")
run_step(
  "retain contigs matching the 40 candidate Ranoidea targets",
  removeOffTargetContigs(
    assembly.directory = reduced_contigs,
    target.markers = target_markers,
    output.directory = target_contigs,
    blast.path = tool_bin,
    memory = memory_gb,
    threads = threads,
    overwrite = overwrite,
    quiet = TRUE
  )
)

count_fasta_records <- function(path) {
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  lines <- readLines(path, warn = FALSE)
  sum(startsWith(lines, ">"))
}

input_summary <- read.delim(
  file.path(example_directory, "provenance", "read_extraction_summary.tsv")
)
input_pairs <- aggregate(
  input_summary$read_pairs_retained,
  by = list(sample = input_summary$sample),
  FUN = sum
)
names(input_pairs)[2] <- "input_pairs"

cleaning_summary <- read.csv(
  file.path(example_directory, "logs", "fastpComplete_summary.csv")
)
cleaned_pairs <- aggregate(
  cleaning_summary$endPairs,
  by = list(sample = cleaning_summary$Sample),
  FUN = sum
)
names(cleaned_pairs)[2] <- "cleaned_pairs"

samples <- sort(unique(input_summary$sample))
output_summary <- merge(input_pairs, cleaned_pairs, by = "sample")
output_summary$draft_contigs <- vapply(
  samples,
  function(sample) count_fasta_records(file.path(draft_contigs, paste0(sample, ".fa"))),
  integer(1)
)
output_summary$nonredundant_contigs <- vapply(
  samples,
  function(sample) count_fasta_records(file.path(reduced_contigs, paste0(sample, ".fa"))),
  integer(1)
)
output_summary$recovered_targets <- vapply(
  samples,
  function(sample) count_fasta_records(file.path(target_contigs, paste0(sample, ".fa"))),
  integer(1)
)
write.csv(
  output_summary,
  file.path(results_directory, "output_summary.csv"),
  row.names = FALSE
)

recovered_target_rows <- do.call(
  rbind,
  lapply(samples, function(sample) {
    path <- file.path(target_contigs, paste0(sample, ".fa"))
    headers <- sub("^>", "", grep("^>", readLines(path, warn = FALSE), value = TRUE))
    data.frame(sample = sample, marker = headers)
  })
)
write.table(
  recovered_target_rows,
  file.path(results_directory, "recovered_targets.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  timings,
  file.path(results_directory, "step_timings.csv"),
  row.names = FALSE
)

environment_lines <- c(
  paste("PhyloProcessR:", as.character(utils::packageVersion("PhyloProcessR"))),
  paste("R:", R.version.string),
  paste("Platform:", R.version$platform),
  paste("Threads:", threads),
  paste("SPAdes memory limit (GB):", memory_gb),
  paste("Tool bin:", tool_bin)
)
writeLines(environment_lines, file.path(results_directory, "run_environment.txt"))

message("\nLaptop example completed. Results: ", results_directory)
