#################################################
## Configuration file for PhyloProcessR
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
working.directory = "/PATH/TO/where/the/stuff/will/happen"
#The processed reads directory
processed.reads = "processed-reads"
#The read folder within processed reads to assemble, pe-merged-reads recommended
assembly.reads = "pe-merged-reads"
# The target markers reference file
target.markers = "PATH/TO/target-markers.fa"

# Global settings
#########################
# number of threads
threads = 16
# Amount of memory to allocate in GB
memory = 120
# TRUE to overwrite previous runs. FALSE the script will resume but will not delete anything.
overwrite = FALSE
# Hide verbose output for each function
quiet = FALSE

# Missing locus recovery settings
#########################
# TRUE = run expandMissingAssembly after the main assembly pipeline
expand.missing = FALSE
# Which read subdirectory within processed.reads to use for Phase 2 mapping.
# Should be paired (non-merged) reads. Options (use whichever is the last step run in workflow 1):
#   "decontaminated-reads" (default — recommended)
#   "cleaned-reads"
#   "error-corrected-reads"
# Do NOT use "pe-merged-reads" — HISAT2 expects paired input
mapping.reads = "decontaminated-reads"
# Reference to use for Phase 2 read mapping:
#   "contig"    = use the best assembled contig from other samples (default; closer match, better read recovery)
#   "reference" = use the original probe/bait sequences (useful when cross-sample contigs are absent or poor)
phase2.reference = "contig"
# TRUE = also attempt to recover loci absent from every sample's assembly,
#        mapping reads directly to the original reference sequences.
#        Can add substantial run time on large datasets.
recover.all.missing = FALSE

#Assembly settings
#########################
#The selected k-mer values for spades
spades.kmer.values = c(33, 55, 77, 99, 127)
#Whether to use mismatch corrector (requires a lot of RAM and resources, recommended if possible)
spades.mismatch.corrector = TRUE
#Whether to save the error-corrected reads produced by SPAdes (ignored if clean.up.spades = TRUE)
save.corrected.reads = FALSE
#TRUE to delete the entire SPAdes working directory for each sample after assembly, keeping only the final .fa file
clean.up.spades = FALSE
# the similarity threshold for redundancy reduction
similarity = 0.95

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/PATH/TO/miniconda3/envs/PhyloProcessR/bin"
cdhit.path    = conda.env
spades.path   = conda.env
blast.path    = conda.env
hisat2.path   = conda.env
samtools.path = conda.env
fastp.path    = conda.env