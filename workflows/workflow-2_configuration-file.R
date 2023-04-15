#################################################
## Configuration file for PhyloProcessR
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory 
working.directory = "/Path/to/where/the/stuff/will/happen"
#The processed reads directory
processed.reads = "processed-reads"
#The read folder within processed reads to assemble, pe-merged-reads recommended
assembly.reads = "pe-merged-reads"

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
#Denovo assembly with spades
denovo.assembly = TRUE
#The selected k-mer values for spades
spades.kmer.values = c(33,55,77,99,127)
#Whether to use mismatch corrector (requires a lot of RAM and resources, recommended if possible)
spades.mismatch.corrector = TRUE
#Whether to save the corrected reads
save.corrected.reads = FALSE
#Clean up large files 
clean.up.spades = FALSE

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/Path/to/conda/env/bin"
fastp.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
spades.path = conda.env
bbnorm.path = conda.env
blast.path = conda.env