source("workflow-1_configuration-file.R")
if (dropbox.download == TRUE && sra.download == TRUE) {
  stop("Workflow 1 accepts one download source at a time. Choose Dropbox or SRA before running.")
}

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
              filter.library.layout   = "PAIRED",
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

# Cleans the reads with one pass of fastp. The fastp command is built from the
# TRUE/FALSE settings in the configuration file, so every step that is TRUE
# runs in that single pass.
if (clean.reads == TRUE) {
  fastpClean(input.reads = input.reads,
             output.directory = paste0(processed.reads, "/cleaned-reads"),
             remove.adaptors = remove.adaptors,
             remove.duplicate.reads = remove.duplicate.reads,
             error.correction = error.correction,
             quality.trim.reads = quality.trim.reads,
             quality.filter = quality.filter,
             low.complexity.filter = low.complexity.filter,
             trim.poly.x = trim.poly.x,
             min.read.length = min.read.length,
             fastp.path = fastp.path,
             threads = threads,
             mem = memory,
             overwrite = overwrite,
             quiet = quiet)
  input.reads = paste0(processed.reads, "/cleaned-reads")
}

# Remove contaminant reads. The function creates references when a list is set.
if (decontamination == TRUE){
  removeContamination(input.reads = input.reads,
                      output.directory = paste0(processed.reads, "/decontaminated-reads"),
                      decontamination.path = if (download.contaminant.genomes) NULL else decontamination.path,
                      decontamination.list = if (download.contaminant.genomes) contaminant.genome.list else NULL,
                      contaminant.directory = "contaminant-references",
                      include.univec = include.univec,
                      include.genbank = contaminant.genbank.accessions,
                      include.fasta = contaminant.fasta,
                      map.match = decontamination.match,
                      samtools.path = samtools.path,
                      bwa.path = bwa.path,
                      threads = threads,
                      mem = memory,
                      overwrite = overwrite,
                      overwrite.reference = overwrite.contaminant.database,
                      overwrite.contaminants = overwrite.contaminant.database,
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
