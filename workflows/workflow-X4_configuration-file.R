#################################################
## Configuration file for workflow X4: novel loci discovery
##
## This workflow discovers novel genomic loci by:
##   1. Identifying reads that do NOT map to any known sequence-capture locus
##   2. Mapping those unmapped reads to a reference genome
##   3. Retaining genomic regions covered by reads in >= min.samples samples
##   4. Assembling per-sample contigs from those regions (SPAdes)
##   5. Filtering for heterozygosity
##   6. Collecting novel contigs (no BLAST — region encoded in contig name)
##   7. Aligning across samples
##   8. Trimming alignments
##
## Prerequisites:
##   - Completed workflow 1 (read preprocessing)
##   - Completed workflow 4 (existing capture alignments in untrimmed_all-markers)
##   - A reference genome FASTA for the organism/group
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
working.directory = "/PATH/TO/PROJECT/DIRECTORY"
# Existing sequence-capture alignments (output of workflow 4 alignment step)
alignment.directory = "data-analysis/alignments/untrimmed_all-markers"
alignment.format = "phylip"
# Directory containing per-sample read subdirectories.
# Point directly to the folder whose immediate children are one directory per sample
# (e.g. the decontaminated-reads subdirectory of your processed reads).
# R1 / R2 files must be present directly inside each sample subdirectory.
read.directory = "processed-reads/decontaminated-reads"
# Full path to reference genome FASTA
genome.file = "/PATH/TO/reference_genome.fa"
# Name for the novel loci dataset (used as prefix for intermediate and output files)
dataset.name = "Novel-Loci"

# Global settings
#########################
# Number of threads
threads = 8
# Amount of memory to allocate in GB
memory = 40
# Whether to overwrite previous runs
overwrite = FALSE
# Hide verbose output for each function
quiet = TRUE

# Novel region discovery settings
#########################
# Minimum number of samples that must cover a region to retain it
min.samples = 4
# Minimum read depth at a site for it to be counted as covered in a sample
min.coverage = 5
# Minimum length in bp for a candidate novel region
min.region.length = 200
# Adjacent covered intervals within this distance (bp) are merged into one region
max.merge.distance = 500
# Minimum MAPQ score — filters reads mapping to repetitive/ambiguous regions
min.mapping.quality = 20

# Assembly settings
#########################
# Minimum reads required per region to attempt SPAdes assembly
min.reads.assemble = 5
# K-mer values for SPAdes
spades.kmer.values = c(33, 55, 77, 99, 127)

# Heterozygosity filter settings
#########################
# TRUE = filter assembled contigs with excessive IUPAC ambiguity before annotation
heterozygote.filter = TRUE
# Maximum allowed proportion of IUPAC ambiguity bases per contig (0-1)
heterozygote.filter.threshold = 0.05
# Contigs shorter than this (bp) are exempt from the proportion filter
heterozygote.min.length = 100

# Contig collection settings
#########################
# TRUE = run collectNovelContigs to consolidate assembled contigs before alignment.
#        No BLAST or probe set is used — the genomic region is already encoded in
#        each contig name by assembleSharedRegions.
annotate.targets = TRUE
# Minimum contig length in bp to retain
min.match.length = 100

# Alignment settings
#########################
# TRUE = run alignTargets after annotation
align.targets = TRUE
# localpair or globalpair alignment algorithm (localpair is slower but more accurate)
alignment.algorithm = "localpair"
# Minimum number of taxa to keep an alignment
min.taxa.alignment = 4
# Maximum pairwise difference from reference sequence for a contig to be removed
removal.threshold = 0.35

# Trimming alignment settings
#########################
# TRUE = run alignment trimming (superTrimmer)
trim.alignments = TRUE
# Minimum alignment length in basepairs to keep an alignment
min.alignment.length = 100
# Run TrimAl to remove highly variable or misaligned columns
run.TrimAl = TRUE
# TRUE = remove samples too divergent from the majority-rule consensus.
# Catches paralogs, off-target captures, and reverse-complemented sequences
# that produce two distinct phylogenetic signals in one alignment.
trim.similarity = TRUE
# Pairwise distance threshold (0-1): samples at or above this distance from
# the consensus are removed. 0.4 removes sequences >40% divergent.
similarity.threshold = 0.4
# TRUE = trim columns below a minimum gap threshold
trim.column = TRUE
# Minimum percent of bases that must be present to keep a column
min.column.gap.percent = 50
# Resolve ambiguous IUPAC sites to a random constituent base
convert.ambiguous.sites = FALSE
# TRUE = trim poorly covered edges of each alignment
trim.external = TRUE
# Minimum percent of bases required to keep an edge column
min.external.percent = 50
# TRUE = remove samples below a minimum coverage threshold
trim.coverage = TRUE
# Minimum percent of bases a sample must have across the alignment
min.coverage.percent = 35
# Minimum number of bases a sample must have across the alignment
min.coverage.bp = 60

# Program paths
#########################
### *** When installing via conda, only the path to the conda bin directory is needed.
### Otherwise modify any of these to the path where the program is found.
conda.env = "PATH/TO/miniconda3/envs/PhyloProcessR/bin"
hisat2.path   = conda.env
samtools.path = conda.env
bedtools.path = conda.env
spades.path   = conda.env
blast.path    = conda.env
mafft.path    = conda.env
trimAl.path   = conda.env

#### End configuration
