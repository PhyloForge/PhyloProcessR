source("workflow-1_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)

##################################################################################################
##################################################################################################
## Step 1: Preprocess reads
##################################################################################################

#Begins by creating processed read directory
dir.create(processed.reads, showWarnings = FALSE)
# The shared results and log directories used by the later workflows
dir.create("data-analysis", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

if (dropbox.download == TRUE){
  #Run download function. The token file is read by the function itself.
  dropboxDownload(sample.spreadsheet = sample.file,
                  dropbox.directory = dropbox.directory,
                  dropbox.token = dropbox.token,
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

# If both sources were active, merge their rename tables into one file so that
# organizeReads sees all samples regardless of origin.
if (dropbox.download == TRUE && sra.download == TRUE) {
  combined.rename = rbind(read.csv("file_rename_dropbox.csv"),
                          read.csv("file_rename_sra.csv"))
  write.csv(combined.rename, "file_rename_combined.csv", row.names = FALSE)
  sample.file = "file_rename_combined.csv"
}

#Organizes reads if scattered elsewhere i.e. creates a sub-dataset
if (organize.reads == TRUE) {
  organizeReads(read.directory = read.directory,
                output.directory = paste0(processed.reads, "/organized-reads"),
                rename.file = sample.file,
                link.reads = link.reads,
                overwrite = overwrite)
  input.reads = paste0(processed.reads, "/organized-reads")
} else {input.reads = read.directory }

if (summary.fastq == TRUE){
  fastqStats(read.directory = input.reads,
             output.name = "logs/fastq-stats",
             read.length = read.length,
             threads = threads,
             mem = memory,
             overwrite = overwrite)
}#end summary.fastq if

# Quick scan of raw reads against the target probe set to flag poor samples early
if (assess.capture == TRUE){
  assessCaptureEfficiency(input.reads = input.reads,
                          output.directory = "data-analysis/sample-capture-assessment",
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

if (remove.adaptors == TRUE && fastp.complete == FALSE) {
  removeAdaptors(input.reads = input.reads,
                 output.directory = paste0(processed.reads, "/adaptor-removed-reads"),
                 fastp.path = fastp.path,
                 threads = threads,
                 mem = memory,
                 overwrite = overwrite,
                 quiet = quiet)
  input.reads = paste0(processed.reads, "/adaptor-removed-reads")
}

# Removes exact PCR duplicates
if (remove.duplicate.reads == TRUE && fastp.complete == FALSE) {
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
if (error.correction == TRUE && fastp.complete == FALSE) {
  readErrorCorrection(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/error-corrected-reads"),
                      fastp.path = fastp.path,
                      threads = threads,
                      mem = memory,
                      overwrite = overwrite,
                      quiet = quiet)
  input.reads = paste0(processed.reads, "/error-corrected-reads")
}

# Trims low quality read ends
if (quality.trim.reads == TRUE && fastp.complete == FALSE) {
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
  #Downloads the contaminant genomes, or uses a local set of genomes
  if (download.contaminant.genomes == TRUE){
    createContaminantDB(decontamination.list = contaminant.genome.list,
                        output.directory = "contaminant-references",
                        include.univec = include.univec,
                        overwrite = overwrite.contaminant.database)
    contaminant.references = "contaminant-references"
  } else {
    if (is.null(decontamination.path) == TRUE){
      stop("Set decontamination.path to a local set of contaminant genomes, or set download.contaminant.genomes = TRUE.")
    }
    contaminant.references = decontamination.path
  }

  ## remove external contamination
  removeContamination(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/decontaminated-reads"),
                      decontamination.path = contaminant.references,
                      map.match = decontamination.match,
                      samtools.path = samtools.path,
                      bwa.path = bwa.path,
                      threads = threads,
                      mem = memory,
                      overwrite = overwrite,
                      overwrite.reference = overwrite.contaminant.database,
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
} #end merge.pe.reads if
