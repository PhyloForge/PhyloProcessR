#################################################
## Configuration file for workflow 4: alignments and trimming
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
# working.directory = "/Path/to/where/the/stuff/is/at"
working.directory = "/Volumes/LaCie/Mantellidae"
# the contig directory to make alignments from
contig.directory = "/Volumes/LaCie/Mantellidae/data-analysis/contigs/5_iupac-contigs"
# The sequence capture target marker file for extraction from contigs
target.file = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
# The name for the dataset
dataset.name = "Descriptive-Name"

#Global settings
#########################
#number of threads
threads = 4
#Amount of memory to allocate in GB
memory = 40
#Whether to overwrite previous runs
overwrite = FALSE
#Resume from previous runs (does not overwrite)
resume = TRUE
#Print verbose output for each function
quiet = TRUE

#Target matching settings
#########################
#Filtering threshold to filter out potential paralogs
heterozyote.filter.threshold = 0.05
#whether to trim to the targets losing flanking sequence (TRUE) or keep the entire contig (FALSE)
trim.to.targets = FALSE
#The minimum match percentage for a contig match to a target
min.match.percent = 60
#The minimum match length in basepairs for a contig match to a target
min.match.length = 40
#The minimum match coverage, contig must overlap by X percent to target
min.match.coverage = 50

# Alignment settings
#########################
# TRUE = to run alignments for the matching targets from above
align.matching.targets = TRUE
# localpair or globalpair, localpair slower but better
alignment.algorithm = "localpair"
# The minimum number of taxa to keep an alignment
min.taxa.alignment = 4
# subset alignments to run multiple instances, uses proportion of targets between 0-1
subset.start = 0
# example, 0, 0.25 to align first quarter, 0.25 to 0.5 to align second
subset.end = 1

# Alignment subsets
#########################
# Trims each alignment to the target marker, leaving out the flanks
trim.to.targets = TRUE
# Trims the target out of each alignment, leaving only the flanks (inverse of previous)
trim.to.flanks = TRUE


# Trimming alignment settings
#########################
#TRUE = to run alignment trimming function batchTrimAlignments
trim.alignments = TRUE
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
#Resolves ambiguous sites to the same arbitrary base
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
### e.g. fastp.path = "/conda/PhyloCap/bin
conda.env = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
blast.path = conda.env
cdhit.path = conda.env
mafft.path = conda.env
trimAl.path = conda.env
julia.path = conda.env

#### End configuration
