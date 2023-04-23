#################################################
## Configuration file for PhyloProcessR
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
#working.directory = "/Path/to/where/the/stuff/will/happen"
working.directory <- "/Volumes/LaCie/Mantellidae"
#The processed reads directory
processed.reads = "processed-reads"
#The read folder within processed reads to assemble, pe-merged-reads recommended
assembly.reads = "pe-merged-reads"
# The target markers reference file
target.markers = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"

# Global settings
#########################
# number of threads
threads = 8
# Amount of memory to allocate in GB
memory = 20
# TRUE to overwrite previous runs. FALSE the script will resume but will not delete anything.
overwrite = FALSE
# Hide verbose output for each function
quiet = FALSE

#Assembly settings
#########################
#The selected k-mer values for spades
spades.kmer.values = c(33,55,77,99,127)
#Whether to use mismatch corrector (requires a lot of RAM and resources, recommended if possible)
spades.mismatch.corrector = TRUE
#Whether to save the corrected reads
save.corrected.reads = FALSE
#Clean up large files 
clean.up.spades = FALSE
# the similarity threshold for redundancy reduction
reduce.redundancy = 0.95

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/Path/to/conda/env/bin"
cdhit.path = conda.env
spades.path = conda.env
blast.path = conda.env