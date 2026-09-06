#################################################
## Configuration file for PhyloProcessR
#################################################

# Package version
#########################
# TRUE installs the latest development/beta version from GitHub before running.
# Keep FALSE for reproducible analyses that use the already installed version.
install.latest.github = FALSE

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
working.directory = "/Path/to/where/the/stuff/will/happen"

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

#Raw read locations
#########################
# The file rename (File, Sample columns) for organizing and setting up names for reads or dropbox download.
# First column is "File" with the file name without read/lane details, second is "Sample" with desired new names. 
# NA will use current file names.
sample.file = "file_rename.csv"
#TRUE = to rename reads and set up for analyses using csv file above. Recommended to keep TRUE. 
organize.reads = TRUE
# TRUE links the organized reads instead of copying them. A link uses no extra
# disk space. Set to FALSE when the reads must be copied, for example to move
# them to a different drive.
link.reads = FALSE
# The input raw read directory, NULL if downloading from dropbox
read.directory = "/Path/to/where/the/raw/reads/are"
#The name for the processed reads folder
processed.reads = "processed-reads"
# TRUE to save a summary csv file of the raw sequence data
summary.fastq = TRUE
# The sequencing read length in base pairs. Used for the megabase pair total.
read.length = 150
# TRUE to map reads to the target probe set and estimate capture efficiency per sample
assess.capture = TRUE
# Path to the target probe/marker FASTA used for sequence capture
target.fasta = "/Path/to/probe-set.fa"

# For downloading reads from dropbox
#########################
#TRUE to download files from personal dropbox folder using token and read path.
dropbox.download = FALSE
#the dropbox directory your files are all contained within if dropbox.download = TRUE
dropbox.directory = "/Dropbox/Path/to/Reads"
# Token file with a saved Dropbox OAuth2 token, saved with saveRDS.
dropbox.token = "/Local/Path/to/Token/token.RDS"
# Skips files not found in file_rename.csv spreadsheet.
skip.not.found = FALSE

# For downloading reads from NCBI SRA (Sequence Read Archive)
#########################
# TRUE to download reads from NCBI SRA using an SraRunInfo CSV file.
# Export your SRA Run Selector results as SraRunInfo.csv from the NCBI SRA
# Run Selector (https://www.ncbi.nlm.nih.gov/Traces/study/) and provide the
# path below. Can be used together with dropbox.download = TRUE. The workflow
# then merges the two rename tables into file_rename_combined.csv.
sra.download = FALSE
# Path to the SraRunInfo CSV downloaded from the NCBI SRA Run Selector.
# Must contain at minimum a 'Run' column with SRR/ERR/DRR accession numbers.
sra.info.file = "SraRunInfo.csv"
# Column in sra.info.file to use directly as sample names.
# NULL (default) auto-builds names from ScientificName_SRRaccession
# (e.g. Hylarana_macrodactyla_SRR11853236). Set to a column name (quoted
# string) to use that column instead.
sra.sample.name.column = NULL
# Only download rows with this LibraryStrategy value (e.g. "Targeted-Capture").
# Set to NULL to download everything in the file regardless of strategy.
sra.filter.strategy = NULL
# Number of retry attempts for a failed file download before skipping.
sra.max.retries = 3
# Seconds to wait between retry attempts.
sra.retry.delay = 10
# TRUE skips samples that fail after all retries (prints a warning).
# FALSE raises an error and stops the run.
sra.skip.not.found = TRUE

#FASTP read cleaning
#########################
# = TRUE to run all processing steps at once (much faster). Overrides settings below.
# = FALSE to run the separate analyses with TRUE below
fastp.complete = TRUE
#TRUE = to run adaptor removal on reads
remove.adaptors = TRUE
#TRUE to remove exact PCR duplicates
remove.duplicate.reads = TRUE
#TRUE to correct errors using the other read pair
error.correction = TRUE

#Other read processing tasks
############################
#Merge paired end reads, helps with assembly
merge.pe.reads = TRUE
# Trims low quality ends off of reads (not recommended, hurts assembly)
quality.trim.reads = FALSE

#Decontamination settings
#########################
#Remove contamination
decontamination = TRUE
# TRUE downloads the contaminant genomes from NCBI into contaminant-references.
# FALSE uses the local genomes in decontamination.path below.
download.contaminant.genomes = TRUE
#The file for the contaminant genomes (Genome, GenBank_Accession columns); only used if download.contaminant.genomes = TRUE
contaminant.genome.list = "decontamination_database.csv"
#A path to a local set of contaminant genomes; only used if download.contaminant.genomes = FALSE
decontamination.path = NULL
#Include the univec contaminant database?
include.univec = TRUE
# Minimum alignment identity (0-1) that makes a read a contaminant; 0.90 = 90% identity.
# A read pair at or above this value is removed and counted. A read pair below
# it is kept.
decontamination.match = 0.90
# TRUE downloads the contaminant genomes again and rebuilds the BWA index. This
# is separate from overwrite so that a new read run does not download every
# contaminant genome again. The index is also rebuilt on its own when the
# contaminant files change.
overwrite.contaminant.database = FALSE

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/Path/to/conda/env/bin"
fastp.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
spades.path = conda.env
blast.path = conda.env
