#################################################
## Configuration file for PhyloProcessR Workflow 3 - Variant Calling
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory where a "dataset.name" directory with variant calling results will be saved.
working.directory = "/Volumes/LaCie/Anax"
# The read directory desired for mapping, recommended "decontaminated-reads". Default shown.
read.directory = "/Volumes/LaCie/Anax/reads"
# The alignment directory desired to have variants called on. Default shown. 
alignment.directory = "/Volumes/LaCie/Anax/data-analysis/alignments/untrimmed_all-markers"
# The name for the dataset
dataset.name = "joint-genotyping"

# Global settings
#########################
# number of threads
threads = 8
# Amount of memory to allocate in GB
memory = 80
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
# TRUE to stop  pipeline when read sets are missing corresponding assemblies; FALSE removes read sets without assemblies
check.assemblies = FALSE
# TRUE to run GATK4 base-recalibrator. Requires high depth. if you observe few SNPs, set use.base.recalibration = FALSE
base.recalibration = FALSE
# TRUE to use the GATK4 base-recalibrator results. Requires high depth. if you observe few SNPs, set this to false
use.base.recalibration = FALSE

# Custom hard filtering thresholds
#########################
# Default GATK4 recommended values are shown here
# For filter explanations see: 
#   https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# Quality score
custom.SNP.QUAL = 30
# Quality by Depth: quality score normalized by depth
custom.SNP.QD = 2
# Strand Odds Ratio: odds ratio of strand bias
custom.SNP.SOR = 3
# Fisher Strand: phred-scaled probability of strand bias
custom.SNP.FS = 60
# Map quality: root mean square mapping quality
custom.SNP.MQ = 40
# Map quality rank sum: compares mapping quality of reads supporting reference and alternative allele
custom.SNP.MQRankSum = -12.5
# Read position rank sum: tests for site position within reads
custom.SNP.ReadPosRankSum = -8
# Indel quality by depth: quality score normalized by depth
custom.INDEL.QD = 2
# Indel quality
custom.INDEL.QUAL = 30
# Indel Fisher strand: phred-scaled probability of strand bias
custom.INDEL.FS = 60
# Indel Read position rank sum: tests for site position within reads
custom.INDEL.ReadPosRankSum = -8

# Output settings
#########################
#VCF file to use. Choices are "SNP", "Indel", "Both". SNP should be used in the majority of cases.
vcf.file = "SNP"
# TRUE to save contigs with ambiguity codes placed at heterozygous sites
ambiguity.codes = TRUE
# TRUE to save contigs using a consensus base (randomly selected) for each heterozygous site.
consensus.sequences = TRUE

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### *** Replace /PATH/TO/ with your system 
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/PATH/TO/miniconda3/envs/PhyloProcessR/bin"
gatk4.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
