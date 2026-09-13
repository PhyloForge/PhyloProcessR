source("workflow-X1_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)

# One flag controls both running base recalibration and using its results. This
# also accepts the older split flags base.recalibration and use.base.recalibration.
old.recalibration = get0("base.recalibration", ifnotfound = NULL)
if (!is.null(old.recalibration)) {
  if (exists("use.base.recalibration") &&
      !identical(as.logical(old.recalibration), as.logical(use.base.recalibration))) {
    stop("base.recalibration and use.base.recalibration disagree. Set only use.base.recalibration.")
  }
  if (!exists("use.base.recalibration")) { use.base.recalibration = old.recalibration }
}
if (!exists("use.base.recalibration")) { use.base.recalibration = FALSE }
if (!exists("ploidy")) { ploidy = 2 }

##################################################################################################
##################################################################################################
## Runs series of functions and organizes results
##################################################################################################

# Begins by creating processed read directory
if (file.exists(paste0("data-analysis/", dataset.name)) == FALSE) {
  dir.create(paste0("data-analysis/", dataset.name))
}#end if

# Shared paths for the dataset. Every stage uses the one dataset-owned reference,
# so all samples share the same contig names, sequences, and coordinates.
mapping.directory = paste0("data-analysis/", dataset.name, "/sample-mapping")
haplotype.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller")
reference.path = paste0("data-analysis/", dataset.name, "/reference/reference.fa")

#Function that prepares the BAM files and sets the metadata correctly for GATK4
prepareBAM(
  read.directory = read.directory,
  output.directory = mapping.directory,
  auto.readgroup = auto.readgroup,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  temp.directory = temp.directory,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

# The expected cohort is every prepared sample. It is captured once and passed to
# every downstream stage so a missing sample fails the run instead of quietly
# reducing the cohort.
expected.samples = list.dirs(mapping.directory, recursive = FALSE, full.names = FALSE)
expected.samples = expected.samples[nzchar(expected.samples)]

# Builds the shared reference and maps every sample against it. reference.mode is
# "consensus" (build from the sample alignments), "target" (use the capture
# target markers file), or "user" (use reference.file).
mapReferenceConsensus(
  mapping.directory = mapping.directory,
  alignment.directory = alignment.directory,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  temp.directory = temp.directory,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet,
  reference.path = reference.path,
  reference.mode = reference.mode,
  target.file = target.file,
  reference.file = reference.file,
  sample.names = expected.samples
)

# Function that calls the haplotypes using GATK4
haplotypeCaller(
  mapping.directory = mapping.directory,
  output.directory = haplotype.directory,
  reference.type = "consensus",
  reference.path = reference.path,
  ploidy = ploidy,
  gatk4.path = gatk4.path,
  temp.directory = temp.directory,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet,
  sample.names = expected.samples
)

# Function that recalibrates bases and calls haplotypes again
if (use.base.recalibration == TRUE) {
  #runs function
  baseRecalibration(
    haplotype.caller.directory = haplotype.directory,
    mapping.directory = mapping.directory,
    gatk4.path = gatk4.path,
    temp.directory = temp.directory,
    threads = threads,
    memory = memory,
    clean.up = clean.up,
    overwrite = overwrite,
    quiet = quiet,
    ploidy = ploidy,
    sample.names = expected.samples,
    reference.path = reference.path
  )

}#end if

# Function that uses GATK4 to genotype and filter samples creating a final VCF of supported SNPs
jointGenotyping(
  haplotype.caller.directory = haplotype.directory,
  output.directory = paste0("data-analysis/", dataset.name, "/genotype-database"),
  use.base.recalibration = use.base.recalibration,
  save.unfiltered = save.unfiltered,
  save.SNPs = save.SNPs,
  save.indels = save.indels,
  save.combined = save.combined,
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
  temp.directory = temp.directory,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet,
  reference.path = reference.path,
  sample.names = expected.samples,
  batch.size = batch.size
)

#END Workflow X1
