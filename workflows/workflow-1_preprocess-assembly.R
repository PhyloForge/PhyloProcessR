#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never", force = TRUE)
library(PhyloCap)
library(foreach)

#setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
#source("/Users/chutter/Dropbox/Research/0_Github/PhyloCap/work-flows/workflow-1_configuration-file.R")
source("workflow-1_configuration-file.R")
setwd(working.directory)

##################################################################################################
##################################################################################################
#################################################
## Step 1: Preprocess reads
##################

##### Ideal processing order
# 1. Organize / Summary
# 2. Remove adaptors (fastp)
# 3. remove duplicates (fastp)
# 4. read correction (fastp)
# 5. trim low quality (fastp)
# 6. Decontamination
# 7. normalize
# 8. Merge PE Reads

#Checks if everything is installed
pass.fail = setupCheck(anaconda.environment =  NULL,
                       fastp.path = fastp.path,
                       samtools.path = samtools.path,
                       bwa.path = bwa.path,
                       spades.path = spades.path,
                       bbnorm.path = bbnorm.path)

if (pass.fail == FALSE){ stop("Some required programs are missing") } else {
  print("all required programs are found, PhyloCap pipeline continuing...")
}

#Begins by creating processed read directory
dir.create("processed-reads")

if (dropbox.download == TRUE){
  #Authorizes token
  rdrop2::drop_auth(rdstoken = dropbox.token)
  #Run download function
  dropboxDownload(sample.spreadsheet = sample.file,
                  dropbox.directory = dropbox.directory,
                  out.directory = paste0(processed.reads, "/raw-reads"),
                  overwrite = overwrite)

  read.directory = paste0(processed.reads, "/raw-reads")
  organize.reads = FALSE
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
  #This function creates a summary of the fastq files per sample for number of reads
  fastqStats(read.directory = input.reads,
             output.name = "fastq-stats",
             read.length = 150,
             threads = threads,
             mem = memory,
             overwrite = overwrite)
}#end summary.fastq if

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
                 resume = resume,
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
                      resume = resume,
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
                      resume = resume,
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
                   resume = resume,
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
                      resume = resume,
                      overwrite = overwrite,
                      overwrite.reference = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/decontaminated-reads")
}


# Normalizes reads
if (normalize.reads == TRUE) {
  normalizeReads(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/normalized-reads"),
                 bbnorm.path = bbnorm.path,
                 threads = threads,
                 mem = memory,
                 resume = resume,
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
                      resume = resume,
                      overwrite = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/pe-merged-reads")
} #end decontamination

dir.create("data-analysis")

if (denovo.assembly == TRUE){
  #Assembles merged paired end reads with spades
  assembleSpades(input.reads = input.reads,
                 output.directory = paste0("data-analysis/spades-assembly"),
                 assembly.directory = "data-analysis/draft-assemblies",
                 mismatch.corrector = spades.mismatch.corrector,
                 kmer.values = spades.kmer.values,
                 threads = threads,
                 memory = memory,
                 overwrite = overwrite,
                 resume = resume,
                 save.corrected.reads = save.corrected.reads,
                 quiet = quiet,
                 spades.path = spades.path)
}#end if

#End script
