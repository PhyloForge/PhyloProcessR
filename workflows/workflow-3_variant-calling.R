#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
source("/Users/chutter/Dropbox/Research/0_Github/R_Projects/PhyloProcessR/workflows/workflow-3_configuration-file.R")
source("workflow-3_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
## Step 1: Preprocess reads
##################################################################################################

#Begins by creating processed read directory
dir.create(dataset.name)

#Function that prepares the BAM files and sets the metadata correctly for GATK4
prepareBAM(
  read.directory = read.directory,
  output.directory = paste0(dataset.name, "/sample-mapping"),
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
mapReference(
  bam.directory = paste0(dataset.name, "/sample-mapping"),
  output.directory = paste0(dataset.name, "/sample-mapping"),
  assembly.directory = assembly.directory,
  check.assemblies = check.assemblies,
  reference.file = reference.file,
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
  bam.directory = paste0(dataset.name, "/sample-mapping"),
  output.directory = paste0(dataset.name, "/haplotype-caller"),
  samtools.path = samtools.path,
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
    haplotype.caller.directory = "variant-discovery/haplotype-caller",
    sample.mapping.directory = "variant-discovery/sample-mapping",
    gatk4.path = gatk4.path,
    threads = threads,
    memory = memory,
    clean.up = clean.up,
    overwrite = overwrite,
    quiet = quiet
  )

}#end if



