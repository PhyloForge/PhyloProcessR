#################################################
## Configuration file for PhyloProcessR Workflow X1 - Joint Genotyping
#################################################

# Package version
#########################
# TRUE installs the latest development/beta version from GitHub before running.
# Keep FALSE for reproducible analyses that use the already installed version.
install.latest.github = FALSE

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory where a "dataset.name" directory with variant calling results will be saved.
working.directory = "/Volumes/LaCie/Anax"
# The read directory desired for mapping, recommended "decontaminated-reads". Default shown.
read.directory = "/Volumes/LaCie/Anax/reads"
# The alignment directory desired to have variants called on. Default shown.
alignment.directory = "/Volumes/LaCie/Anax/data-analysis/alignments/untrimmed_all-markers"

# Shared reference source. All samples map and are genotyped against one reference.
#   "consensus" = one majority consensus per alignment in alignment.directory (default).
#   "target"    = use the capture target markers file (target.file).
#   "user"      = use a reference FASTA you supply (reference.file).
reference.mode = "consensus"
# The capture target markers FASTA, used only when reference.mode = "target".
target.file = "/PATH/TO/marker-seqs.fa"
# A user-supplied reference FASTA, used only when reference.mode = "user".
reference.file = "/PATH/TO/reference.fa"
# Temporary directory where temporary files are saved
temp.directory = working.directory
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
# TRUE runs GATK4 base recalibration (BQSR) and uses its recalibrated GVCFs. One
# flag controls both steps. Requires high depth; keep FALSE if you observe few SNPs.
use.base.recalibration = FALSE
# Ploidy passed to both the initial and the recalibrated haplotype caller.
ploidy = 2
# Number of samples imported per GenomicsDB batch. Bounds import memory only;
# every sample still enters every locus.
batch.size = 50
# TRUE to save unfiltered variant calling data (recommended)
save.unfiltered = TRUE
# TRUE to save SNPs vcf separately
save.SNPs = TRUE
# TRUE to save indels vcf separately
save.indels = TRUE
# TRUE to save a combined SNP and indels vcf separately
save.combined = TRUE

# Custom hard filtering thresholds
#########################
# Default GATK4 recommended values are shown here. These are cohort-level record
# FILTER expressions on site annotations. A passing record does not guarantee that
# every sample genotype at that site has adequate depth; this workflow applies no
# per-sample genotype-depth filter and writes variant-only VCFs.
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

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### *** Replace /PATH/TO/ with your system
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
gatk4.path = conda.env
samtools.path = conda.env
bwa.path = conda.env
