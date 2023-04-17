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
  reference.file = reference.file,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)







if (dropbox.download == TRUE){
  #Authorizes token
  rdrop2::drop_auth(rdstoken = dropbox.token)
  #Run download function
  dropboxDownload(sample.spreadsheet = sample.file,
                  dropbox.directory = dropbox.directory,
                  output.directory = paste0(processed.reads, "/raw-reads"),
                  overwrite = overwrite,
                  skip.not.found = skip.not.found)

  read.directory = paste0(processed.reads, "/raw-reads")
  organize.reads = TRUE
  sample.file = "file_rename_dropbox.csv"
}#end if

#Organizes reads if scattered elsewhere i.e. creates a sub-dataset
if (organize.reads == TRUE) {
  organizeReads(read.directory = read.directory,
                output.dir = paste0(processed.reads, "/organized-reads"),
                rename.file = sample.file,
                overwrite = overwrite)
  input.reads = paste0(processed.reads, "/organized-reads")
} else {input.reads = read.directory }

# if (summary.fastq == TRUE){
#   #This function creates a summary of the fastq files per sample for number of reads
#   fastqStats(read.directory = input.reads,
#              output.name = "fastq-stats",
#              read.length = 150,
#              threads = threads,
#              mem = memory,
#              overwrite = overwrite)
# }#end summary.fastq if

#The complete processing through fastp at once. +++ for speed.
if (fastp.complete == TRUE) {
  fastpComplete(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/cleaned-reads"),
                 fastp.path = fastp.path,
                 threads = threads,
                 mem = memory,
                 overwrite = overwrite,
                 quiet = quiet)
  input.reads = paste0(processed.reads, "/cleaned-reads")
}

if (remove.adaptors == TRUE & fastp.complete == FALSE) {
  removeAdaptors(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/adaptor-removed-reads"),
                 fastp.path = fastp.path,
                 threads = threads,
                 mem = memory,
                 overwrite = overwrite,
                 quiet = quiet)
  input.reads = paste0(processed.reads, "/adaptor-removed-reads")
}

# Runs read error correction
if (remove.duplicate.reads == TRUE & fastp.complete == FALSE) {
  removeDuplicateReads(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/deduped-reads"),
                      fastp.path = fastp.path,
                      threads = threads,
                      mem = memory,
                      overwrite = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/deduped-reads")
}

# Runs read error correction
if (error.correction == TRUE & fastp.complete == FALSE) {
  readErrorCorrection(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/error-corrected-reads"),
                      fastp.path = fastp.path,
                      threads = threads,
                      mem = memory,
                      overwrite = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/error-corrected-reads")
}

# Normalizes reads
if (quality.trim.reads == TRUE & fastp.complete == FALSE) {
  qualityTrimReads(input.reads = input.reads,
                   output.directory = paste0(processed.reads, "/quality-trimmed-reads"),
                   fastp.path = fastp.path,
                   threads = threads,
                   mem = memory,
                   overwrite = overwrite,
                   quiet = quiet)
  input.reads = paste0(processed.reads, "/quality-trimmed-reads")
}

#Runs decontamination of reads
if (decontamination == TRUE){
  #Creates the database by downloading
  createContaminantDB(decontamination.list = contaminant.genome.list,
                      output.directory = "contaminant-references",
                      include.human = include.human,
                      include.univec = include.univec,
                      overwrite = overwrite)

  ## remove external contamination
  removeContamination(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/decontaminated-reads"),
                      decontamination.path = "contaminant-references",
                      map.match = decontamination.match,
                      samtools.path = samtools.path,
                      bwa.path = bwa.path,
                      threads = threads,
                      mem = memory,
                      overwrite = overwrite,
                      overwrite.reference = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/decontaminated-reads")
}

# Normalizes reads
if (normalize.reads == TRUE) {
  #Normalizes reads using ORNA
  normalizeReads(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/normalized-reads"),
                 orna.path = orna.path,
                 threads = threads,
                 memory = memory,
                 overwrite = overwrite,
                 quiet = quiet)
  input.reads = paste0(processed.reads, "/normalized-reads")
}

#merge paired-end reads
if (merge.pe.reads == TRUE){
  #merge paired end reads
  mergePairedEndReads(input.reads = input.reads,
                      output.directory =  paste0(processed.reads, "/pe-merged-reads"),
                      fastp.path = fastp.path,
                      threads = threads,
                      mem = memory,
                      overwrite = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/pe-merged-reads")
} #end decontamination
