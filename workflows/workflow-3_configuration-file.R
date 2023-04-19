#################################################
## Configuration file for PhyloProcessR Workflow 3 - Variant Calling
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory where a "dataset.name" directory with variant calling results will be saved.
working.directory = "/Volumes/LaCie/Mantellidae"
# The read directory desired for mapping, recommended "decontaminated-reads"
read.directory = "/Volumes/LaCie/Mantellidae/reads"
# The assembly directory desired to have variants called on. 
assembly.directory = "/Volumes/LaCie/Mantellidae/expanded-assemblies"
#The name for the dataset
dataset.name = "variant-calling"

# Global settings
#########################
# number of threads
threads = 6
# Amount of memory to allocate in GB
memory = 20
# TRUE to overwrite previous runs. FALSE the script will resume but will not delete anything.
overwrite = FALSE
# Hide verbose output for each function
quiet = FALSE
# deletes intermediate files
clean.up = TRUE

# Variant calling pipeline settings
#########################
# TRUE to determine and name read groups from Illumina headers. FALSE to give arbitrary names.
auto.readgroup = TRUE
# TRUE to stop the pipeline when read sets do not have corresponding assemblies; FALSE removes read sets without assemblies
check.assemblies = FALSE
# TRUE to run the GATK4 base-recalibrator. 
# Uses initial strongly supported haplotype calls, refines base calls
# Recommended only if you have sufficient sequencing depth greater than 20X. 
base.recalibration = TRUE

# Output settings
#########################
#VCF file to use. Choices are "SNP", "Indel", "Both". SNP should be used in the majority of cases.
vcf.file = "SNP"
# TRUE to save contigs with ambiguity codes placed at heterozygous sites
ambuiguity.codes = TRUE
# TRUE to save contigs using a consensus base (randomly selected) for each heterozygous site.
consensus.sequences = TRUE

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
gatk4.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
