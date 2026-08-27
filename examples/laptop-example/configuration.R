# Laptop-scale PhyloProcessR example configuration.
#
# Optional environment variables:
#   PHYLOPROCESSR_BIN        Conda bin directory containing external programs.
#   PHYLOPROCESSR_THREADS    CPU threads (default: 2).
#   PHYLOPROCESSR_MEMORY_GB  SPAdes memory limit in GB (default: 4).
#   PHYLOPROCESSR_OVERWRITE  TRUE/FALSE (default: FALSE).

threads <- as.integer(Sys.getenv("PHYLOPROCESSR_THREADS", unset = "2"))
memory_gb <- as.integer(Sys.getenv("PHYLOPROCESSR_MEMORY_GB", unset = "4"))
overwrite <- tolower(Sys.getenv("PHYLOPROCESSR_OVERWRITE", unset = "false")) %in%
  c("true", "t", "1", "yes")

tool_bin <- Sys.getenv("PHYLOPROCESSR_BIN", unset = "")
if (!nzchar(tool_bin)) {
  detected_fastp <- Sys.which("fastp")
  if (nzchar(detected_fastp)) {
    tool_bin <- dirname(detected_fastp)
  }
}

if (!nzchar(tool_bin) || !dir.exists(tool_bin)) {
  stop(
    "Set PHYLOPROCESSR_BIN to the bin directory of the PhyloProcessR ",
    "Conda environment (normally $CONDA_PREFIX/bin)."
  )
}

raw_reads <- file.path(example_directory, "raw-reads")
target_markers <- file.path(
  example_directory,
  "targets",
  "Ranoidea-V1_forty-candidate-locus.fa"
)
results_directory <- file.path(example_directory, "results")

spades_kmers <- c(21, 33, 55)
spades_mismatch_corrector <- FALSE
reduce_redundancy_similarity <- 0.95
