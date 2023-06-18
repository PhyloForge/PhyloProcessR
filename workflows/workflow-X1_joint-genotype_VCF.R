#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

source("/Users/chutter/Dropbox/Research/0_Github/R_Projects/PhyloProcessR/workflows/workflow-X1_configuration-file.R")
source("workflow-X1_configuration-file.R")

setwd(working.directory)

##################################################################################################
##################################################################################################
## Runs series of functions and organizes results
##################################################################################################

# Begins by creating processed read directory
if (file.exists(paste0("data-analysis/", dataset.name)) == FALSE) {
  dir.create(paste0("data-analysis/", dataset.name))
}#end if

#Function that prepares the BAM files and sets the metadata correctly for GATK4
prepareBAM(
  read.directory = read.directory,
  output.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  auto.readgroup = auto.readgroup,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

#Function that prepares the BAM files and sets the metadata correctly for GATK4
mapReferenceConsensus(
  mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  alignment.directory = alignment.directory,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

# Function that calls the haplotypes using GATK4
haplotypeCaller(
  mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  output.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller"),
  reference.type = "consensus",
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

# Function that recalibrates bases and calls haplotypes again
if (base.recalibration == TRUE) {
  #runs function
  baseRecalibration(
    haplotype.caller.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller"),
    mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
    gatk4.path = gatk4.path,
    threads = threads,
    memory = memory,
    clean.up = clean.up,
    overwrite = overwrite,
    quiet = quiet
  )

}#end if

# Function that uses GATK4 to genotype and filter samples creating a final VCF of supported SNPs
genotypeSamples(
  mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  haplotype.caller.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller"),
  output.directory = paste0("data-analysis/", dataset.name, "/sample-genotypes"),
  use.base.recalibration = use.base.recalibration,
  custom.SNP.QD =  custom.SNP.QD,
  custom.SNP.QUAL =  custom.SNP.QUAL,
  custom.SNP.SOR =  custom.SNP.SOR,
  custom.SNP.FS =  custom.SNP.FS,
  custom.SNP.MQ =  custom.SNP.MQ,
  custom.SNP.MQRankSum =  custom.SNP.MQRankSum,
  custom.SNP.ReadPosRankSum =  custom.SNP.ReadPosRankSum,
  custom.INDEL.QD =  custom.INDEL.QD,
  custom.INDEL.QUAL =  custom.INDEL.QUAL,
  custom.INDEL.FS =  custom.INDEL.FS,
  custom.INDEL.ReadPosRankSum =  custom.INDEL.ReadPosRankSum,  
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

if (consensus.sequences == TRUE) {
  # Function that converts SNP files back into finished and SNP called contigs, choose format
  VCFtoContigs(
    genotype.directory = paste0("data-analysis/", dataset.name, "/sample-genotypes"),
    mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
    output.directory = paste0("data-analysis/contigs/4_consensus-contigs"),
    vcf.file = vcf.file,
    consensus.sequences = TRUE,
    ambiguity.codes = FALSE,
    gatk4.path = gatk4.path,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet
  )
}

if (ambiguity.codes == TRUE) {
  # Function that converts SNP files back into finished and SNP called contigs, choose format
  VCFtoContigs(
    genotype.directory = paste0("data-analysis/", dataset.name, "/sample-genotypes"),
    mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
    output.directory = paste0("data-analysis/contigs/5_iupac-contigs"),
    vcf.file = vcf.file,
    consensus.sequences = FALSE,
    ambiguity.codes = TRUE,
    gatk4.path = gatk4.path,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet
  )
}

#END Workflow 3