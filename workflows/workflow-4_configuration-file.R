#################################################
## Configuration file for workflow 4: alignments and trimming
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
working.directory = "/PATH/TO/PROJECT/DIRECTORY"
# the contig directory to make alignments from
contig.directory = "/data-analysis/contigs/5_iupac-contigs"
# The sequence capture target marker file for extraction from contigs
target.file = "marker-seqs.fa"
# The name for the dataset
dataset.name = "Descriptive-Name"

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

#Target annotation settings
#########################
# TRUE to run target annotation, FALSE if its already been done to skip
annotate.targets = TRUE
#Filtering threshold to filter out potential paralogs
heterozyote.filter.threshold = 0.05
#The minimum match percentage for a contig match to a target
min.match.percent = 60
#The minimum match length in basepairs for a contig match to a target
min.match.length = 40
#The minimum match coverage, contig must overlap by X percent to target
min.match.coverage = 50

# Alignment settings
#########################
# TRUE = to run alignments for the annotated targets from above
align.targets = TRUE
# localpair or globalpair, localpair slower but better
alignment.algorithm = "localpair"
# The minimum number of taxa to keep an alignment
min.taxa.alignment = 4
# The pairwise difference threshold from reference for removal
removal.threshold = 0.35
# subset alignments to run multiple instances, uses proportion of targets between 0-1
subset.start = 0
# example, 0, 0.25 to align first quarter, 0.25 to 0.5 to align second
subset.end = 1

#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
conda.env = "PATH/TO/miniconda3/envs/PhyloProcessR/bin"
blast.path = conda.env
cdhit.path = conda.env
mafft.path = conda.env

#### End configuration
