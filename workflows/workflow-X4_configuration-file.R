#################################################
## Configuration file for workflow X4: novel loci discovery
##
## This workflow discovers novel genomic loci by:
##   1. Identifying reads that do NOT map to any known sequence-capture locus
##   2. Mapping those unmapped reads to a reference genome
##   3. Retaining genomic regions covered by reads in >= min.samples samples
##   4. Assembling per-sample contigs from those regions (SPAdes)
##   5. Filtering for heterozygosity
##   6. Annotating and aligning across samples (as in workflow 4)
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
# Processed reads directory (contains per-sample subdirectories)
processed.reads = "processed-reads"
# Read subdirectory to use for mapping — must be paired (non-merged)
mapping.reads = "decontaminated-reads"
# Full path to reference genome FASTA
genome.file = "/PATH/TO/reference_genome.fa"
# Name for the novel loci dataset (used for intermediate files in annotateTargets)
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

# Annotation settings
#########################
# TRUE = run annotateTargets to match contigs to novel target sequences
annotate.targets = TRUE
# Minimum BLAST match percent identity
min.match.percent = 60
# Minimum BLAST match length in bp
min.match.length = 40
# Minimum contig overlap with the target (percent)
min.match.coverage = 50
# TRUE = retain rather than remove potential paralogs
retain.paralogs = FALSE

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

# Program paths
#########################
### *** When installing via conda, only the path to the conda bin directory is needed.
### Otherwise modify any of these to the path where the program is found.
conda.env = "PATH/TO/miniconda3/envs/PhyloProcessR/bin"
hisat2.path    = conda.env
samtools.path  = conda.env
bedtools.path  = conda.env
spades.path    = conda.env
blast.path     = conda.env
cdhit.path     = conda.env
mafft.path     = conda.env

#### End configuration
