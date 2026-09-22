source("workflow-3_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)

variant.directory = file.path("data-analysis", dataset.name)
contig.directory = file.path("data-analysis", "contigs")
depth.directory = file.path("data-analysis", "depth")

##################################################################################################
##################################################################################################
## Runs series of functions and organizes results
##################################################################################################

# Create the directory for variant-calling intermediate files.
if (dir.exists(variant.directory) == FALSE) {
  dir.create(variant.directory, recursive = TRUE)
}#end if

# The GATK temporary directory must not be a general data directory because the
# workflow removes it after a successful run.
if (length(temp.directory) != 1 || is.na(temp.directory) ||
    nchar(temp.directory) == 0) {
  stop("temp.directory must name one dedicated temporary directory.")
}
temp.path = normalizePath(temp.directory, mustWork = FALSE)
working.path = normalizePath(working.directory, mustWork = TRUE)
unsafe.temp.paths = c(normalizePath("/", mustWork = TRUE),
                      normalizePath(path.expand("~"), mustWork = TRUE),
                      working.path)
if (temp.path %in% unsafe.temp.paths) {
  stop("temp.directory must not be the filesystem root, home directory, or working.directory.")
}

# The dedicated GATK temporary directory must exist before any GATK call.
dir.create(temp.directory, showWarnings = FALSE, recursive = TRUE)

#Function that prepares the BAM files and sets the metadata correctly for GATK4
prepareBAM(
  read.directory = read.directory,
  output.directory = file.path(variant.directory, "sample-mapping"),
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

#Function that maps each sample to its own assembly
retained.samples = mapReferenceSample(
  mapping.directory = file.path(variant.directory, "sample-mapping"),
  assembly.directory = assembly.directory,
  check.assemblies = check.assemblies,
  samtools.path = samtools.path,
  bwa.path = bwa.path,
  gatk4.path = gatk4.path,
  temp.directory = temp.directory,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet
)

# Function that calls the haplotypes using GATK4
haplotypeCaller(
  mapping.directory = file.path(variant.directory, "sample-mapping"),
  output.directory = file.path(variant.directory, "haplotype-caller"),
  reference.type = "sample",
  ploidy = ploidy,
  gatk4.path = gatk4.path,
  temp.directory = temp.directory,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet,
  sample.names = retained.samples
)

# Function that recalibrates bases and calls haplotypes again.
# Only runs when use.base.recalibration = TRUE; the results feed directly into
# genotypeSamples below via the same flag.
if (use.base.recalibration == TRUE) {
  baseRecalibration(
    haplotype.caller.directory = file.path(variant.directory, "haplotype-caller"),
    mapping.directory = file.path(variant.directory, "sample-mapping"),
    gatk4.path = gatk4.path,
    temp.directory = temp.directory,
    threads = threads,
    memory = memory,
    clean.up = clean.up,
    overwrite = overwrite,
    quiet = quiet,
    ploidy = ploidy,
    sample.names = retained.samples
  )
}#end if

# Function that uses GATK4 to genotype and filter samples creating a final VCF of supported SNPs
genotypeSamples(
  mapping.directory = file.path(variant.directory, "sample-mapping"),
  haplotype.caller.directory = file.path(variant.directory, "haplotype-caller"),
  output.directory = file.path(variant.directory, "sample-genotypes"),
  use.base.recalibration = use.base.recalibration,
  temp.directory = temp.directory,
  custom.SNP.QD = custom.SNP.QD,
  custom.SNP.QUAL = custom.SNP.QUAL,
  custom.SNP.SOR = custom.SNP.SOR,
  custom.SNP.FS = custom.SNP.FS,
  custom.SNP.MQ = custom.SNP.MQ,
  custom.SNP.MQRankSum = custom.SNP.MQRankSum,
  custom.SNP.ReadPosRankSum = custom.SNP.ReadPosRankSum,
  custom.INDEL.QD = custom.INDEL.QD,
  custom.INDEL.QUAL = custom.INDEL.QUAL,
  custom.INDEL.FS = custom.INDEL.FS,
  custom.INDEL.ReadPosRankSum = custom.INDEL.ReadPosRankSum,
  gatk4.path = gatk4.path,
  threads = threads,
  memory = memory,
  overwrite = overwrite,
  quiet = quiet,
  sample.names = retained.samples
)

depth.files = NULL
if (depth.filter.mode != "none" || !is.null(max.n.proportion)) {
  depth.files = calculateSampleDepth(
    mapping.directory = file.path(variant.directory, "sample-mapping"),
    output.directory = depth.directory,
    sample.names = retained.samples,
    use.base.recalibration = use.base.recalibration,
    samtools.path = samtools.path,
    overwrite = overwrite,
    quiet = quiet
  )
}

if (consensus.sequences == TRUE) {
  # Function that converts SNP files back into finished and SNP called contigs, choose format
  VCFtoContigs(
    genotype.directory = file.path(variant.directory, "sample-genotypes"),
    mapping.directory = file.path(variant.directory, "sample-mapping"),
    output.directory = file.path(contig.directory, "4_consensus-contigs"),
    vcf.file = vcf.file,
    consensus.sequences = TRUE,
    ambiguity.codes = FALSE,
    temp.directory = temp.directory,
    gatk4.path = gatk4.path,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet,
    sample.names = retained.samples,
    depth.files = depth.files,
    depth.filter.mode = depth.filter.mode,
    min.site.depth = min.site.depth,
    min.mean.depth = min.mean.depth,
    max.n.proportion = max.n.proportion,
    use.base.recalibration = use.base.recalibration,
    samtools.path = samtools.path,
    ploidy = ploidy
  )
}

if (ambiguity.codes == TRUE) {
  # Function that converts SNP files back into finished and SNP called contigs, choose format
  VCFtoContigs(
    genotype.directory = file.path(variant.directory, "sample-genotypes"),
    mapping.directory = file.path(variant.directory, "sample-mapping"),
    output.directory = file.path(contig.directory, "5_iupac-contigs"),
    vcf.file = vcf.file,
    consensus.sequences = FALSE,
    ambiguity.codes = TRUE,
    temp.directory = temp.directory,
    gatk4.path = gatk4.path,
    threads = threads,
    memory = memory,
    overwrite = overwrite,
    quiet = quiet,
    sample.names = retained.samples,
    depth.files = depth.files,
    depth.filter.mode = depth.filter.mode,
    min.site.depth = min.site.depth,
    min.mean.depth = min.mean.depth,
    max.n.proportion = max.n.proportion,
    use.base.recalibration = use.base.recalibration,
    samtools.path = samtools.path,
    ploidy = ploidy
  )
}

# Keep temporary files after a failure for diagnosis. Remove them only after
# every requested workflow step finishes successfully.
if (dir.exists(temp.directory) == TRUE) {
  unlink(temp.directory, recursive = TRUE, force = TRUE)
}
if (dir.exists(temp.directory) == TRUE) {
  warning("Workflow 3 completed, but could not remove temp.directory: ",
          temp.directory)
}

#END Workflow 3
