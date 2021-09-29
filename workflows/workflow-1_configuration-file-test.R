#################################################
## Configuration file for PhyloCap
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
#working.directory = "/home/c111h652/scratch/Centrolenidae"
working.directory = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset"
#The file for the contaminant genomes (Genome, Accession columns); use NULL if download.contaminant.genomes = F
contaminant.genome.list = "decontamination_database.csv"
#The name for the dataset
dataset.name = "Test"

#Raw read locations
#########################
#The file rename (File, Sample columns) for organizing reads or dropbox download
sample.file = NA
#TRUE = to subset reads from a larger pool using the "sample.file" above
organize.reads = FALSE
#The input raw read directory
read.directory = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/raw-reads/"
#The name for the processed reads folder
processed.reads = "processed-reads"
#TRUE to download files from personal dropbox folder using token and read path
dropbox.download = FALSE
#the dropbox directory your files are all contained within
dropbox.directory = "/Research/3_Sequence-Database/Raw-Reads"
#dropbox.tok = "/Users/chutter/Dropbox/dropbox-token.RDS"
dropbox.token = "/home/c111h652/dropbox-token.RDS"
#TRUE to save a summary csv file of the raw sequence data
summary.fastq = TRUE

#Global settings
#########################
#number of threads
threads = 16
#Amount of memory to allocate in GB
memory = 120
#Whether to overwrite previous runs
overwrite = FALSE
#Resume from previous runs (does not overwrite)
resume = TRUE
#Print verbose output for each function
quiet = TRUE

#FASTP read cleaning
#########################
# = TRUE to run all processing steps at once (much faster). Overrides settings below.
fastp.complete = TRUE
#TRUE = to run adaptor removal on reads
remove.adaptors = TRUE
#TRUE to remove exact PCR duplicates
remove.duplicate.reads = TRUE
#TRUE to correct errors using the other read pair
error.correction = TRUE
#Trims low quality ends off of reads
quality.trim.reads = TRUE
#normalizes read depth to facilitate assembly
normalize.reads = TRUE
#Merge paired end reads
merge.pe.reads = TRUE

#Decontamination settings
#########################
#Remove contamination
decontamination = TRUE
#Download contaminat genomes from genbank if TRUE
download.contaminant.genomes = TRUE
#Alternatively, a path can be set to a downloaded set of contaminant genomes if downloading does not work; NULL if downloading
decontamination.path = NULL
#include the human genome? Unless human is study organism or UCEs in mammals are used
include.human = TRUE
#Include the univec contaminant database?
include.univec = TRUE
#The matching proportion of bases for a contaminant hit and removal
decontamination.match = 0.99

#Assembly settings
#########################
#Denovo assembly with spades
denovo.assembly = TRUE
#The selected k-mer values for spades
spades.kmer.values = c(33,55,77,99,127)
#Whether to use mismatch corrector (requires a lot of RAM and resources)
spades.mismatch.corrector = TRUE
#Whether to save the corrected reads
save.corrected.reads = FALSE
#Clean up large files *** ADD INTO ANALYSES
clean.up.spades = FALSE


#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
### e.g. fastp.path = "/conda/PhyloCap/bin
#conda.env = "/panfs/pfs.local/work/bi/c111h652/conda/PhyloCap/bin/"
conda.env = "/Users/chutter/miniconda3/bin/"
fastp.path = conda.env
samtools.path = conda.env
bwa.path = "/usr/local/bin"
spades.path = conda.env
bbnorm.path = "/usr/local/bin"
blast.path = conda.env

