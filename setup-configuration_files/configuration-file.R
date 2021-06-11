#################################################
## Configuration file for PhyloCap
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
work.dir = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset"
#The file rename (File, Sample columns) for organizing reads. Set to NULL if not needed
file.rename = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/file_rename.csv"
#The file for the contaminant genomes (Genome, Accession columns); use NULL if download.contaminant.genomes = F
contaminant.genome.list = "/Users/chutter/Dropbox/Research/0_Github/PhyloCap/setup-configuration_files/decontamination_database.csv"
#The sequence capture target marker file for extraction from contigs
target.file = "/Users/chutter/Dropbox/Research/0_Github/PhyloCap/setup-configuration_files/Ranoidea_All-Markers_Apr21-2019.fa"
#The input raw read directory
read.dir = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/raw-reads"
#The name for the dataset
dataset.name = "Test"
#The name for the processed reads folder
processed.reads = "processed-reads"

#Global settings
#########################
#number of threads
threads = 4
#Amount of memory to allocate in GB
memory = 8
#Whether to overwrite previous runs
overwrite = FALSE
#Resume from previous runs (does not overwrite)
resume = TRUE
#Print verbose output for each function
quiet = TRUE

#Pre-processing settings
#########################
#TRUE = to organize reads and rename to sample names from the "file_rename.csv" above
organize.reads = TRUE
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
spades.kmer.values = c(21,33,55,77,99,127)
#Whether to use mismatch corrector (requires a lot of RAM and resources)
spades.mismatch.corrector = TRUE
#Whether to save the corrected reads
save.corrected.reads = FALSE
#Clean up large files *** ADD INTO ANALYSES
clean.up.spades = FALSE


#Target matching settings
#########################
#TRUE = to run target matching on contigs
match.targets = TRUE
#Directory of assembled contigs to be used for target matching
assembly.directory = "draft-assemblies"
#whether to trim to the targets losing flanking sequence (TRUE) or keep the entire contig (FALSE)
trim.to.targets = FALSE
#The directory of contigs to use for target matching
#contig.directory = "draft-scaffolds" #TO ADD IN
#The minimum match percentage for a contig match to a target
min.match.percent = 60
#The minimum match length for a contig match to a target
min.match.length = 40
#The minimum match coverage, contig must overlap by X proportion to target
min.match.coverage = 50

#Alignment settings
#########################
#TRUE = to run alignments for the matching targets from above
align.matching.targets = TRUE
#The minimum number of taxa to keep an alignment
min.taxa.alignment = 4

#Alignment trimming settings
#########################
#TRUE = to run alignment trimming function batchTrimAlignments
trim.alignments = TRUE
#The minimum number of taxa to keep an alignment
min.taxa.alignment.trim = 5
#The minimum alignment basepairs to keep an alignment
min.alignment.length = 100
#The maximum gaps from throughout the entire alignment to keep an alignment
max.alignment.gap.percent = 50
#run the trimming program TAPER to trim out unalignment sample segments
run.TAPER = TRUE
#run the trimming program TrimAl to remove high variable or misaligned columns
run.TrimAl = TRUE
#Whether to trim out columns below a certain threshold
trim.column = TRUE
#The percent of bases that must be present to keep a column
min.column.gap.percent = 50
#Resolves ambiguous sites to the same arbitrary base
convert.ambiguous.sites = TRUE
#TRUE = to output an alignment assessment spreadsheet
alignment.assess = TRUE
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

#Tree settings
#########################
#TRUE = to estimate gene trees for each alignment
estimate.gene.trees = TRUE
#The minimum number of taxa to keep a gene tree
min.taxa.tree = 4
#TRUE = Clean up iqtree files except the ML tree
cleanup.genetrees = TRUE #Only saves the ML tree; FALSE saves all IQTree files for each gene tree

#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
### e.g. fastp.path = "/conda/PhyloCap/bin
fastp.path = "/Users/chutter/conda/PhyloCap/bin"
samtools.path = "/Users/chutter/conda/PhyloCap/bin"
bwa.path = "/Users/chutter/conda/PhyloCap/bin"
spades.path = "/Users/chutter/conda/PhyloCap/bin"
bbmap.path = "/Users/chutter/conda/PhyloCap/bin"
blast.path = "/Users/chutter/conda/PhyloCap/bin"
mafft.path = "/Users/chutter/conda/PhyloCap/bin"
iqtree.path = "/Users/chutter/conda/PhyloCap/bin"
trimAl.path = "/Users/chutter/conda/PhyloCap/bin"
taper.path = "/Users/chutter/conda/PhyloCap/bin"
julia.path = "/Users/chutter/conda/PhyloCap/bin"

#### End configuration
