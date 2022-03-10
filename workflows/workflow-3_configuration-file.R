################################################
## Configuration file for PhyloCap
#################################################

#Directories and input files
#########################
# *** Full paths should be used whenever possible
#The main working directory
#working.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Guibemantis-liber"
#The file for the contaminant genomes (Genome, Accession columns); use NULL if download.contaminant.genomes = F
#alignment.directory = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Guibemantis-liber/alignments/trimmed_all-markers"

#Cluster
working.directory = "/home/c111h652/scratch/Microhylidae"
#The file for the contaminant genomes (Genome, Accession columns); use NULL if download.contaminant.genomes = F
alignment.directory = "/home/c111h652/scratch/Microhylidae/data-analysis/alignments/trimmed_all-markers"
#The name for the alignment dataset
dataset.name = "trimmed_all-markers"

#Global settings
#########################
#number of threads
threads = 16
#Amount of memory to allocate in GB
memory = 120
#Whether to overwrite previous runs
overwrite = FALSE
#Print verbose output for each function
quiet = TRUE

#Gene tree estimation settings
#########################
#The minimum number of taxa to keep an alignment
min.taxa.tree = 5
#Select between "file" (by gene from concatenating) or "merge" which is automatic partition finding
partition.scheme = "merge"
#Whether to partition by codon position (only if partition.scheme = "file")
codon.partition = FALSE
#Program to use (only IQTree available for now )
program = "IQTREE"
#nuclear or mitochondrial data
msub.type = "nuclear"
#number of ultra-fast bootstrap replicates, must be greater than 1000
uf.bootstrap = 1000
#number of top models to assess for sequence evolution models. Lower number = faster
rcluster = 10
#clean up gene trees, deletes intermediate files except ML best tree
cleanup.genetrees = TRUE

#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
### e.g. fastp.path = "/conda/PhyloCap/bin
conda.env = "/panfs/pfs.local/work/bi/c111h652/conda/PhyloCap/bin/"
iqtree.path = conda.env

#### End configuration
