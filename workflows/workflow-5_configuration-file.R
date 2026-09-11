#################################################
## Configuration file for workflow 5: trimming
#################################################

# Package version
#########################
# TRUE installs the latest development/beta version from GitHub before running.
# Keep FALSE for reproducible analyses that use the already installed version.
install.latest.github = FALSE

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
working.directory = "/PATH/TO/PROJECT/DIRECTORY"
# The sequence capture target marker file for extraction from contigs
target.file = "marker-seqs.fa"
# Gene metadata file. Required columns are "marker" and "gene".
# "Marker" and "Gene" are also accepted for compatibility.
feature.gene.names = "data-analysis/gene_metadata.txt"

# Global settings
#########################
#number of threads
threads = 8
#Amount of memory to allocate in GB
memory = 40
#Whether to overwrite previous runs
overwrite = FALSE
#Print verbose output for each function
quiet = TRUE

# MACSE Exon alignment refinement
#########################
# TRUE = run MACSE on all no-flank alignments and save separate coding outputs.
# Use this only when all input targets are coding and in the correct reading frame.
# These outputs do not replace the standard no-flank unlinked dataset.
run.macse = TRUE
# The genetic code to use for MACSE (default: 1 for standard nuclear, 2 for vertebrate mitochondrial)
macse.genetic.code = 1

# Alignment subset
#########################
# TRUE = run makeAlignmentSubset on subset.alignment.directory
run.subset = FALSE
# Existing alignment directory to subset. The default is produced only when
# concatenate.genes = FALSE and trim.alignments = TRUE.
subset.alignment.directory = "data-analysis/alignments/trimmed_all-markers"
# Name used for the output subdirectory: data-analysis/alignments/<subset.name>
subset.name = "subset_markers"
# Path to a fasta file whose sequence names identify the alignments to keep
subset.fasta = NULL
# Regular expression used when subset.reference = "grep".
subset.grep.string = NULL
# Method used to match alignments: "fasta" (name matching), "grep" (pattern), or "blast" (similarity)
subset.reference = "fasta"

# Novel markers integration
#########################
# TRUE = incorporate novel loci from workflow X4 (untrimmed_novel-markers) into
#        trimming and unlinked datasets. Expects that directory to exist.
include.novel.markers = FALSE

# Alignment subsets
#########################
# Concatenates exons from same gene
concatenate.genes = TRUE
# minimum number of exons needed to make a concatenated gene
minimum.exons = 2
# Gathers all unlinked alignments i.e. concatenated genes and single exon genes and UCEs
gather.unlinked = TRUE
# Trims each alignment to the target marker, leaving out the flanks
trim.to.targets = TRUE
# Trims the target out of each alignment, leaving only the flanks (inverse of previous)
trim.to.flanks = TRUE

# Trimming alignment settings
#########################
# TRUE = run superTrimmer on the full-marker dataset.
# Target-only and flank-only dataset construction use their own switches.
trim.alignments = TRUE
# The minimum number of taxa. An alignment at this value is rejected.
min.taxa.alignment = 4
# The minimum alignment length. Final assessment rejects an alignment at this value.
min.alignment.length = 100
#The maximum gaps from throughout the entire alignment to keep an alignment
max.alignment.gap.percent = 50
#run the trimming program TrimAl to remove high variable or misaligned columns
run.TrimAl = TRUE
# TRUE = remove samples too divergent from the majority-rule consensus.
# Catches paralogs, off-target captures, and reverse-complemented sequences
# that produce two distinct phylogenetic signals in one alignment.
trim.similarity = TRUE
# Pairwise distance threshold (0-1): samples at or above this distance from
# the consensus are removed. 0.4 removes sequences >40% divergent.
similarity.threshold = 0.4
#Whether to trim out columns below a certain threshold
trim.column = TRUE
# Gap percentage at which a column is removed. For example, 30 removes columns
# with 30 percent gaps or more.
min.column.gap.percent = 50
# Resolves ambiguous IUPAC sites with the deterministic A/T-priority mapping.
convert.ambiguous.sites = FALSE
#TRUE = to externally trim alignment edges
trim.external = TRUE
#The minimum percent of bases that must be present to keep a column on the edges
min.external.percent = 50
# TRUE = remove samples below the coverage thresholds. Percentage coverage is
# relative to the longest sample, not the full alignment width.
trim.coverage = TRUE
#The minimum percent of bases that must be present to keep a sample
min.coverage.percent = 35
#The minimum number of bases that must be present to keep a sample
min.coverage.bp = 60
# TRUE = to output an alignment assessment spreadsheet and filter alignments
alignment.assess = TRUE

#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
conda.env = "PATH/TO/miniconda3/envs/PhyloProcessR/bin"
blast.path = conda.env
mafft.path = conda.env
trimAl.path = conda.env
macse.path = conda.env

#### End configuration
