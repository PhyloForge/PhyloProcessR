source("workflow-X2_configuration-file.R")
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
## Workflow X2: capture efficiency assessment on cleaned reads
## Run workflow 1 first to produce the cleaned reads. This workflow points the
## read-statistic and capture-assessment functions at that read directory.
##################################################################################################

# Confirm the cleaned reads from workflow 1 exist before doing any work
if (dir.exists(read.directory) == FALSE ||
    length(list.files(read.directory, recursive = TRUE)) == 0) {
  stop("No reads found in read.directory: ", read.directory,
       ". Run workflow 1 first to generate the cleaned reads needed for workflow X2.")
}

dir.create("logs", showWarnings = FALSE)

# Quick per-sample FastQ statistics on the cleaned reads
if (summary.fastq == TRUE) {
  fastqStats(read.directory = read.directory,
             output.name = "logs/fastq-stats",
             read.length = read.length,
             threads = threads,
             mem = memory,
             overwrite = overwrite)
}

# Maps cleaned reads against the target probe set to estimate capture efficiency
assessCaptureEfficiency(input.reads = read.directory,
                        output.directory = "data-analysis/sample-capture-assessment",
                        target.fasta = target.fasta,
                        bwa.path = bwa.path,
                        samtools.path = samtools.path,
                        threads = threads,
                        mem = memory,
                        overwrite = overwrite,
                        quiet = quiet)

# Optional barcode identification on the cleaned reads
if (run.barcode.scan == TRUE) {
  if (!requireNamespace("MItoTrawlR", quietly = TRUE)) {
    stop("run.barcode.scan = TRUE requires the MItoTrawlR package.")
  }
  MItoTrawlR::barcodeSampleScan(input.reads = read.directory,
                                output.directory = "barcode-assessment",
                                barcode.fasta = barcode.fasta,
                                database.fasta = barcode.database.fasta,
                                hits.per.sample = 5,
                                per.max.length = barcode.per.max.length,
                                min.iterations = barcode.min.iterations,
                                max.iterations = barcode.max.iterations,
                                min.ref.id = barcode.min.ref.id,
                                bbmap.path = bbmap.path,
                                spades.path = spades.path,
                                cap3.path = cap3.path,
                                blast.path = blast.path,
                                memory = memory,
                                threads = threads,
                                overwrite = overwrite,
                                quiet = quiet)
}

### End workflow X2
