#################################################
## Configuration file for PhyloProcessR
## Workflow X2: capture efficiency assessment on cleaned reads
## Run workflow 1 first to produce the cleaned reads. This workflow points the
## read-statistic and capture-assessment functions at that read directory.
#################################################

# Package version
#########################
# TRUE installs the latest development/beta version from GitHub before running.
# Keep FALSE for reproducible analyses that use the already installed version.
install.latest.github = FALSE

# Working directory
#########################
# *** Full paths should be used whenever possible
working.directory = "/Path/to/where/the/stuff/will/happen"

# Global settings
#########################
# Number of threads
threads = 16
# Amount of memory to allocate in GB
memory = 120
# Hide verbose output for each function
quiet = FALSE
# Whether to overwrite previous runs
overwrite = FALSE

# Read input settings
#########################
# Name of the folder for all processed reads (created by workflow 1)
processed.reads = "processed-reads"
# Directory of cleaned reads to assess. Defaults to the cleaned-reads output of
# workflow 1. Point it elsewhere if the reads are in another location.
read.directory = paste0(processed.reads, "/cleaned-reads")
# Expected read length in bp, used for MegaBasePairs calculation
read.length = 150
# TRUE to run per-sample FastQ statistics on the reads
summary.fastq = TRUE

# Capture assessment settings
#########################
# Full path to the target probe/marker FASTA used for sequence capture
target.fasta = "/Path/to/probe-set.fa"

# Barcode scan settings  (calls MItoTrawlR::barcodeSampleScan)
#########################
# TRUE = run barcode identification on the reads
# FALSE = skip barcode identification entirely
run.barcode.scan = TRUE
# Full path to a FASTA of barcode reference sequence(s) used to recruit reads
# (e.g. a 16S rRNA or COI representative for the target taxa)
barcode.fasta = "/Path/to/barcode-reference.fa"
# Full path to a FASTA of named barcode sequences for local BLAST identification.
# NULL (recommended) queries NCBI nt remotely — no local database needed,
# but requires an internet connection. Remote queries are slower; for large
# datasets set this to a curated local database to avoid NCBI rate limits.
barcode.database.fasta = NULL
# Minimum number of iterative assembly rounds before convergence is tested.
# 3 is appropriate for barcode regions (much shorter than a mitogenome).
barcode.min.iterations = 3
# Maximum number of iterative assembly rounds.
barcode.max.iterations = 10
# Starting BBMap minimum identity for read recruitment. Permissive on first pass
# (0.70) to cast a wide net; automatically tightened to 0.95 once contigs assemble.
barcode.min.ref.id = 0.70
# Fraction above reference length that triggers the max-length guard.
barcode.per.max.length = 0.50

# Program paths
#########################
# When installed via conda, only the path to the conda bin directory is needed
conda.env = "/Path/to/conda/env/bin"
samtools.path = conda.env
bwa.path = conda.env
bbmap.path = conda.env
spades.path = conda.env
cap3.path = conda.env
blast.path = conda.env
