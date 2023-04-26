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

#Assembly settings
#########################
#The selected k-mer values for spades
spades.kmer.values = c(33, 55, 77, 99, 127)
#Whether to use mismatch corrector (requires a lot of RAM and resources, recommended if possible)
spades.mismatch.corrector = TRUE
#Whether to save the corrected reads
save.corrected.reads = FALSE
#Clean up large files 
clean.up.spades = FALSE
# the similarity threshold for redundancy reduction
similarity = 0.95

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/PATH/TO/miniconda3/envs/PhyloProcessR/bin"
cdhit.path = conda.env
spades.path = conda.env
blast.path = conda.env