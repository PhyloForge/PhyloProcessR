source("workflow-X0_configuration-file.R")

if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }

  remotes::install_github("PhyloForge/PhyloProcessR",
                          upgrade = "never",
                          dependencies = FALSE)
}

library(PhyloProcessR)
setwd(working.directory)

##################################################################################################
##################################################################################################
## Workflow X0: Rapid per-sample read screening
##################################################################################################

screenReads(
  read.directory = read.directory,
  processed.reads = processed.reads,
  use.dropbox = use.dropbox,
  sample.file = sample.file,
  dropbox.directory = dropbox.directory,
  dropbox.token = dropbox.token,
  delete.raw.reads = delete.raw.reads,
  delete.cleaned.reads = delete.cleaned.reads,
  target.fasta = target.fasta,
  read.length = read.length,
  remove.adaptors = fastp.remove.adaptors,
  remove.duplicate.reads = fastp.remove.duplicate.reads,
  error.correction = fastp.error.correction,
  quality.trim.reads = fastp.quality.trim.reads,
  quality.filter = fastp.quality.filter,
  low.complexity.filter = fastp.low.complexity.filter,
  trim.poly.x = fastp.trim.poly.x,
  min.read.length = fastp.min.read.length,
  run.barcode.scan = run.barcode.scan,
  barcode.fasta = barcode.fasta,
  barcode.database.fasta = barcode.database.fasta,
  barcode.hits.per.sample = barcode.hits.per.sample,
  barcode.min.iterations = barcode.min.iterations,
  barcode.max.iterations = barcode.max.iterations,
  barcode.min.ref.id = barcode.min.ref.id,
  barcode.per.max.length = barcode.per.max.length,
  fastp.path = fastp.path,
  bwa.path = bwa.path,
  samtools.path = samtools.path,
  bbmap.path = bbmap.path,
  spades.path = spades.path,
  cap3.path = cap3.path,
  blast.path = blast.path,
  threads = threads,
  memory = memory,
  quiet = quiet
)

### End workflow X0
