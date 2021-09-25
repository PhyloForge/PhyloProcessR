#################################################
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
working.directory = "/home/c111h652/scratch/Guibemantis"
#The file for the contaminant genomes (Genome, Accession columns); use NULL if download.contaminant.genomes = F
alignment.directory = "/home/c111h652/scratch/Guibemantis/data-analysis/alignments/trimmed_all-markers"
#The name for the dataset
dataset.name = "guibemantis"

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

#Filter settings
#########################
filter.length = NULL
filter.sample = c(0.5, 0.7)
filter.prop.pis = NULL
filter.count.pis = NULL

#Concatenation settings
#########################
#select between "concatenated" and "directory"
#concatenated: a single phylip file of the concatenated loci and also a partition file
#directory: the filtered set of markers will be saved to a separate folder
output.type = "concatenated"
#The minimum number of taxa to keep an alignment
min.alignments = 5
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
#number of top models to assess for sequence evolution models
rcluster = 10

#Program paths
#########################
### *** Modify any of these from NULL to the path that the program is found if R is not detecting system paths
### e.g. fastp.path = "/conda/PhyloCap/bin
conda.env = "/panfs/pfs.local/work/bi/c111h652/conda/PhyloCap/bin/"
#conda.env = "/usr/local/bin/iqtree2"
iqtree.path = conda.env

#### End configuration
