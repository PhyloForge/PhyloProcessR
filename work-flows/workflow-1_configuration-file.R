#################################################
## Configuration file for PhyloCap
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
working.directory = "/home/c111h652/scratch/Centrolenidae"
#The file rename (File, Sample columns) for organizing reads. Set to NULL if not needed
sample.file = "Centrolenidae.csv"
#The file for the contaminant genomes (Genome, Accession columns); use NULL if download.contaminant.genomes = F
contaminant.genome.list = "decontamination_database.csv"
#The sequence capture target marker file for extraction from contigs
target.file = "Hyloidea_All-Markers_Apr21-2019.fa"
#The name for the dataset
dataset.name = "Centrolenidae"

#Raw read locations
#########################
#The input raw read directory
read.dir = NULL
#The name for the processed reads folder
processed.reads = "processed-reads"
#TRUE to download files from personal dropbox folder using token and read path
dropbox.download = TRUE
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

#Pre-processing settings
#########################
#TRUE = to organize reads and rename to sample names from the "file_rename.csv" above
organize.reads = FALSE
#TRUE = to run adaptor removal on reads
remove.adaptors = TRUE
#Merge paired end reads
merge.pe.reads = TRUE
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
conda.env = "/panfs/pfs.local/work/bi/c111h652/conda/PhyloCap/bin/"
fastp.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
spades.path = conda.env
bbmap.path = conda.env
blast.path = conda.env

