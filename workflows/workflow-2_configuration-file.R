#################################################
## Configuration file for workflow -2: alignments and trimming
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
work.dir = "/home/c111h652/scratch/Mantidactylus"
#The file rename (File, Sample columns) for organizing reads. Set to NULL if not needed
sample.file = "file_rename.csv"
#The file for the contaminant genomes (Genome, Accession columns); use NULL if download.contaminant.genomes = F
contaminant.genome.list = "decontamination_database.csv"
#The sequence capture target marker file for extraction from contigs
target.file = "Ranoidea_All-Markers_Apr21-2019.fa"
#The input raw read directory
read.dir = "/home/c111h652/scratch/All_Frogs/novogene/raw_data"
#The name for the dataset
dataset.name = "Mantidactylus"
#The name for the processed reads folder
processed.reads = "processed-reads"

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

#Target matching settings
#########################
#TRUE = to run target matching on contigs
match.targets = TRUE
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
min.taxa.alignment.trim = 4
#The minimum alignment basepairs to keep an alignment
min.alignment.length = 100
#The maximum gaps from throughout the entire alignment to keep an alignment
max.alignment.gap.percent = 50
#run the trimming program TAPER to trim out unalignment sample segments
run.TAPER = FALSE
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

#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
### e.g. fastp.path = "/conda/PhyloCap/bin
conda.env = "/panfs/pfs.local/work/bi/c111h652/conda/PhyloCap/bin/"
bbmap.path = conda.env
blast.path = conda.env
mafft.path = conda.env
trimAl.path = conda.env
taper.path = conda.env
julia.path = "~/programs/julia/bin/julia"

#### End configuration
