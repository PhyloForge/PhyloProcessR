#################################################
## Configuration file for PhyloProcessR Workflow 3 - Variant Calling
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory, must already exist. "data-analysis" from the assembly step recommended. 
working.directory = "/Volumes/LaCie/Mantellidae/data-analysis"
# The read directory desired for mapping, recommended "decontaminated-reads"
read.directory = "/Volumes/LaCie/Mantellidae/reads"
# The assembly directory desired to have variants called on. 
assembly.directory = "/Volumes/LaCie/Mantellidae/expanded-assemblies"
# The assembly directory desired to have variants called on.
reference.file = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
#The name for the dataset
dataset.name = "variant-calling"

# Global settings
#########################
# number of threads
threads = 6
# Amount of memory to allocate in GB
memory = 20
# TRUE to overwrite previous runs. FALSE the script will resume but will not delete anything.
overwrite = TRUE
# Hide verbose output for each function
quiet = FALSE

# BAM and read mapping settings
#########################
#TRUE to determine and name read groups from Illumina headers. FALSE to give arbitrary names. 
auto.readgroup = TRUE



  # Debugging
  # library(PhyloCap)
  # library(foreach)
  # setwd("/Volumes/LaCie/Mantellidae")
  # assembly.directory <- "/Volumes/LaCie/Mantellidae/expanded-assemblies"
  # output.directory <- "variant-discovery/sample-mapping"
  # reference.file <- "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
  # bam.directory <- "/Volumes/LaCie/Mantellidae/variant-discovery/sample-mapping"

  # iterations <- 5
  # gatk4.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # bwa.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # hisat2.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"

  # auto.readgroup <- T
  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE





#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
gatk4.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
blast.path = conda.env