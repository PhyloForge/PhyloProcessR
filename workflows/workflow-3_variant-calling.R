#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
source("/Users/chutter/Dropbox/Research/0_Github/R_Projects/PhyloProcessR/workflows/workflow-3_configuration-file.R")
source("workflow-3_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
## Runs series of functions and organizes results
##################################################################################################

# Begins by creating processed read directory
dir.create(paste0("data-analysis/", dataset.name))

#Function that prepares the BAM files and sets the metadata correctly for GATK4
prepareBAM(
  read.directory = read.directory,
  assembly.directory = assembly.directory,
  output.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  auto.readgroup = auto.readgroup,
  check.assemblies = check.assemblies,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

#Function that prepares the BAM files and sets the metadata correctly for GATK4
mapReference(
  mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  assembly.directory = assembly.directory,
  check.assemblies = check.assemblies,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

# Function that calls the haplotypes using GATK4
haplotypeCallerGATK4(
  mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  output.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller"),
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

# Function that recalibrates bases and calls haplotypes again
if (base.recalibration == TRUE){
  #runs function
  baseRecalibration(
    haplotype.caller.directory = paste0(dataset.name, "/haplotype-caller"),
    mapping.directory = paste0(dataset.name, "/sample-mapping"),
    gatk4.path = gatk4.path,
    threads = threads,
    memory = memory,
    clean.up = clean.up,
    overwrite = overwrite,
    quiet = quiet
  )

}#end if

# Function that uses GATK4 to genotype and filter samples creating a final VCF of supported SNPs
variants.genotypeSamples(
  mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping"),
  haplotype.caller.directory = paste0(dataset.name, "/haplotype-caller"),
  output.directory = paste0("data-analysis/", dataset.name, "/sample-genotypes"),
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
    output.directory = paste0("data-analysis/", dataset.name, "/iupac-contigs"),
    vcf.file = vcf.file,
    consensus.sequences = TRUE,
    ambiguity.codes = FALSE,
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
    output.directory = paste0("data-analysis/", dataset.name, "/consensus-contigs"),
    vcf.file = vcf.file,
    consensus.sequences = FALSE,
    ambiguity.codes = TRUE,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet
  )
}

#END Workflow 3