#################################################
## Configuration file for workflow 5: trimming
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
working.directory = "/PATH/TO/PROJECT/DIRECTORY"
# The sequence capture target marker file for extraction from contigs
target.file = "marker-seqs.fa"
#feature gene name metadata file, column one: "Marker"; column: "Gene"
feature.gene.names = "/data-analysis/gene_metadata.txt"

#Global settings
#########################
#number of threads
threads = 8
#Amount of memory to allocate in GB
memory = 40
#Whether to overwrite previous runs
overwrite = FALSE
#Print verbose output for each function
quiet = TRUE

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
# TRUE = to run alignment trimming function batchTrimAlignments
trim.alignments = TRUE
# The minimum number of taxa to keep an alignment
min.taxa.alignment = 4
#The minimum alignment basepairs to keep an alignment
min.alignment.length = 100
#The maximum gaps from throughout the entire alignment to keep an alignment
max.alignment.gap.percent = 50
#run the trimming program TrimAl to remove high variable or misaligned columns
run.TrimAl = TRUE
#Whether to trim out columns below a certain threshold
trim.column = TRUE
#The percent of bases that must be present to keep a column
min.column.gap.percent = 50
#Resolves ambiguous sites to random base it could be
convert.ambiguous.sites = FALSE
#TRUE = to externally trim alignment edges
trim.external = TRUE
#The minimum percent of bases that must be present to keep a column on the edges
min.external.percent = 50
#TRUE = to trim samples below a certain coverage (percent bases present out of total alignment) threshold
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

#### End configuration
