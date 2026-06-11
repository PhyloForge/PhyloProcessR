#Installs updated package version
devtools::install_github("chutter/PhyloProcessR", upgrade = "never", dependencies = FALSE)
library(PhyloProcessR)

#source("/Users/chutter/Dropbox/Research/0_Github/PhyloProcessR/work-flows/workflow-1_configuration-file.R")
source("workflow-1_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
## Step 1: Preprocess reads
##################################################################################################

#Begins by creating processed read directory
dir.create(processed.reads)

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

# Download reads directly from NCBI SRA via ENA HTTPS mirrors.
# Provide sra.info.file = path to SraRunInfo.csv from the NCBI SRA Run Selector.
if (sra.download == TRUE){
  sraDownload(sra.info.file          = sra.info.file,
              sample.name.column      = sra.sample.name.column,
              output.directory        = paste0(processed.reads, "/raw-reads"),
              filter.library.strategy = sra.filter.strategy,
              max.retries             = sra.max.retries,
              retry.delay             = sra.retry.delay,
              skip.not.found          = sra.skip.not.found,
              overwrite               = overwrite,
              quiet                   = quiet)

  read.directory = paste0(processed.reads, "/raw-reads")
  organize.reads = TRUE
  sample.file    = "file_rename_sra.csv"
}#end if

#Organizes reads if scattered elsewhere i.e. creates a sub-dataset
if (organize.reads == TRUE) {
  organizeReads(read.directory = read.directory,
                output.dir = paste0(processed.reads, "/organized-reads"),
                rename.file = sample.file,
                overwrite = overwrite)
  input.reads = paste0(processed.reads, "/organized-reads")
} else {input.reads = read.directory }

if (summary.fastq == TRUE){
  fastqStats(read.directory = input.reads,
             output.name = "fastq-stats",
             read.length = 150,
             threads = threads,
             mem = memory,
             overwrite = overwrite)
}#end summary.fastq if

# Quick scan of raw reads against the target probe set to flag poor samples early
if (assess.capture == TRUE){
  assessCaptureEfficiency(input.reads = input.reads,
                          output.directory = "sample-capture-assessment",
                          target.fasta = target.fasta,
                          bwa.path = bwa.path,
                          samtools.path = samtools.path,
                          threads = threads,
                          mem = memory,
                          overwrite = overwrite,
                          quiet = quiet)
}#end assess.capture if

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
