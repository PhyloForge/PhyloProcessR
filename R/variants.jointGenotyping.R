#' @title jointGenotyping
#'
#' @description Performs joint genotyping across a fixed cohort of samples at each
#'   locus using GATK4. Every expected sample must have exactly one GVCF (the
#'   BQSR or the standard caller output, depending on use.base.recalibration) and
#'   its index; a missing sample stops the run rather than silently genotyping a
#'   subset. For each locus the cohort GVCFs are combined into an owned GenomicsDB
#'   workspace, GenotypeGVCFs genotypes the combined data, SNPs and indels are
#'   separated, hard-filtered, merged, and written to per-locus VCF files. Loci
#'   are processed in parallel and every external command is checked. A cohort
#'   record ties the results to the exact sample set, reference, BQSR selection,
#'   and filter thresholds; a changed cohort or reference requires overwrite.
#'
#' @param haplotype.caller.directory path to the directory of per-sample GVCF
#'   sub-directories (output of haplotypeCaller()).
#' @param output.directory path where the GenomicsDB workspace and per-locus VCF
#'   subdirectories are created.
#' @param use.base.recalibration logical; if TRUE the BQSR-recalibrated GVCFs are
#'   used as input.
#' @param save.unfiltered logical; if FALSE unfiltered VCF subdirectories are
#'   deleted after filtering.
#' @param save.SNPs logical; if FALSE SNP-specific VCF subdirectories are deleted.
#' @param save.indels logical; if FALSE indel-specific VCF subdirectories are
#'   deleted.
#' @param save.combined logical; if FALSE the combined (all variants) VCF
#'   subdirectories are deleted.
#' @param custom.SNP.QD,custom.SNP.QUAL,custom.SNP.SOR,custom.SNP.FS,custom.SNP.MQ,custom.SNP.MQRankSum,custom.SNP.ReadPosRankSum
#'   numeric hard-filter thresholds for SNPs.
#' @param custom.INDEL.QD,custom.INDEL.QUAL,custom.INDEL.FS,custom.INDEL.ReadPosRankSum
#'   numeric hard-filter thresholds for indels.
#' @param gatk4.path GATK executable directory/path or NULL to search the PATH.
#' @param temp.directory GATK JVM temp directory; NULL uses a temporary directory.
#' @param threads number of loci to process in parallel.
#' @param memory total JVM heap budget in GB across concurrent loci.
#' @param overwrite logical; if TRUE the output directory is rebuilt.
#' @param quiet logical; if TRUE tool output is suppressed while logs are retained.
#' @param reference.path shared reference FASTA used for genotyping and filtering.
#' @param sample.names expected cohort sample IDs; NULL uses every GVCF
#'   sub-directory.
#' @param batch.size number of samples imported per GenomicsDB batch. Every
#'   sample still enters every locus; this only bounds memory during import.
#'
#' @details The custom SNP and indel thresholds are cohort-level record FILTER
#'   expressions (QUAL, QD, MQ, and related site annotations). A record that
#'   passes these filters does not guarantee that every sample genotype at that
#'   site has adequate depth. This function applies no per-sample genotype-depth
#'   filter and produces variant-only VCFs, not per-sample FASTA or N-masked
#'   sequences. Absence of a record is therefore not a measure of per-sample
#'   coverage or callability.
#'
#' @return invisibly the loci names; writes per-locus VCF files to
#'   output.directory subdirectories.
#'
#' @export

jointGenotyping = function(haplotype.caller.directory = "haplotype-caller",
                          output.directory = "genotype-database",
                          use.base.recalibration = FALSE,
                          save.unfiltered = TRUE,
                          save.SNPs = TRUE,
                          save.indels = TRUE,
                          save.combined = TRUE,
                          custom.SNP.QD = 2,
                          custom.SNP.QUAL = 30,
                          custom.SNP.SOR = 3,
                          custom.SNP.FS = 60,
                          custom.SNP.MQ = 40,
                          custom.SNP.MQRankSum = -12.5,
                          custom.SNP.ReadPosRankSum = -8,
                          custom.INDEL.QD = 2,
                          custom.INDEL.QUAL = 30,
                          custom.INDEL.FS = 60,
                          custom.INDEL.ReadPosRankSum = -8,
                          gatk4.path = NULL,
                          temp.directory = NULL,
                          threads = 1,
                          memory = 1,
                          overwrite = FALSE,
                          quiet = TRUE,
                          reference.path = "index/reference.fa",
                          sample.names = NULL,
                          batch.size = 50) {

  # Quick checks
  if (is.null(haplotype.caller.directory) || !dir.exists(haplotype.caller.directory)) {
    stop("Haplotype caller directory not found.")
  }
  if (!file.exists(reference.path)) { stop("Reference not found: ", reference.path) }
  reference.dict = sub("\\.fa$", ".dict", reference.path)
  if (!all(file.exists(c(paste0(reference.path, ".fai"), reference.dict)))) {
    stop("Reference index (.fai) and dictionary (.dict) are required next to ", reference.path)
  }

  # At least one final variant product must be requested
  if (save.SNPs == FALSE && save.indels == FALSE && save.combined == FALSE) {
    stop("At least one of save.SNPs, save.indels, or save.combined must be TRUE.")
  }

  # Filter thresholds must be single finite numbers. Negative rank-sum values are
  # allowed.
  thresholds = list(custom.SNP.QD = custom.SNP.QD, custom.SNP.QUAL = custom.SNP.QUAL,
                    custom.SNP.SOR = custom.SNP.SOR, custom.SNP.FS = custom.SNP.FS,
                    custom.SNP.MQ = custom.SNP.MQ, custom.SNP.MQRankSum = custom.SNP.MQRankSum,
                    custom.SNP.ReadPosRankSum = custom.SNP.ReadPosRankSum,
                    custom.INDEL.QD = custom.INDEL.QD, custom.INDEL.QUAL = custom.INDEL.QUAL,
                    custom.INDEL.FS = custom.INDEL.FS, custom.INDEL.ReadPosRankSum = custom.INDEL.ReadPosRankSum)
  for (name in names(thresholds)) {
    value = thresholds[[name]]
    if (length(value) != 1 || !is.finite(value)) { stop(name, " must be a single finite number.") }
  }
  if (length(batch.size) != 1 || !is.finite(batch.size) || batch.size < 1 ||
      batch.size != as.integer(batch.size)) {
    stop("batch.size must be a positive integer.")
  }

  if (is.null(temp.directory)) { temp.directory = tempdir() }
  .ensureDirectory(temp.directory, "temporary directory")
  .ensureDirectory("logs/sample_logs", "sample log directory")

  # Resolves exactly one GVCF and index for every expected sample
  gvcf.name = if (use.base.recalibration) "gatk4-bqsr-haplotype-caller.g.vcf.gz" else "gatk4-haplotype-caller.g.vcf.gz"
  discovered = list.dirs(haplotype.caller.directory, recursive = FALSE, full.names = FALSE)
  discovered = discovered[nzchar(discovered)]
  if (is.null(sample.names)) { sample.names = discovered }
  if (length(sample.names) == 0) { stop("The cohort is empty; no samples to genotype.") }

  gvcf.files = file.path(haplotype.caller.directory, sample.names, gvcf.name)
  gvcf.index = paste0(gvcf.files, ".tbi")
  missing = sample.names[!file.exists(gvcf.files) | !file.exists(gvcf.index)]
  if (length(missing) > 0) {
    stop("Missing ", gvcf.name, " or its index for sample(s): ", paste(missing, collapse = ", "))
  }

  # Confirms one unique sample per GVCF and compatible contigs against the reference
  reference.seq = Biostrings::readDNAStringSet(reference.path)
  loci.names = names(reference.seq)
  vcf.sample.ids = character(length(sample.names))
  for (i in seq_along(sample.names)) {
    header = .vcfHeaderLines(gvcf.files[i])
    ids = .vcfSamples(header)
    if (length(ids) != 1) {
      stop("Expected one sample in ", gvcf.files[i], " but found ", length(ids), ".")
    }
    contigs = .vcfContigs(header)
    if (length(contigs) > 0 && !all(loci.names %in% contigs)) {
      stop("The GVCF for ", sample.names[i], " is missing reference loci; it may use a different reference.")
    }
    vcf.sample.ids[i] = ids
  }
  if (anyDuplicated(vcf.sample.ids)) {
    stop("Duplicate sample IDs across the cohort GVCFs: ",
         paste(unique(vcf.sample.ids[duplicated(vcf.sample.ids)]), collapse = ", "))
  }

  # Describes the cohort so a changed request cannot mix versions across loci
  filter.settings = list(SNP.QD = custom.SNP.QD, SNP.QUAL = custom.SNP.QUAL,
                         SNP.SOR = custom.SNP.SOR, SNP.FS = custom.SNP.FS, SNP.MQ = custom.SNP.MQ,
                         SNP.MQRankSum = custom.SNP.MQRankSum, SNP.ReadPosRankSum = custom.SNP.ReadPosRankSum,
                         INDEL.QD = custom.INDEL.QD, INDEL.QUAL = custom.INDEL.QUAL,
                         INDEL.FS = custom.INDEL.FS, INDEL.ReadPosRankSum = custom.INDEL.ReadPosRankSum)
  info = file.info(gvcf.files)
  current.cohort = list(samples = sort(vcf.sample.ids),
                        reference = tools::md5sum(reference.path),
                        use.base.recalibration = use.base.recalibration,
                        gvcf.size = info$size[order(vcf.sample.ids)],
                        gvcf.mtime = as.character(info$mtime[order(vcf.sample.ids)]),
                        save = c(save.unfiltered, save.SNPs, save.indels, save.combined),
                        filters = filter.settings)
  cohort.path = file.path(output.directory, "cohort-record.rds")

  if (overwrite == TRUE && dir.exists(output.directory)) {
    unlink(output.directory, recursive = TRUE)
  }
  if (file.exists(cohort.path) && overwrite == FALSE) {
    if (!identical(readRDS(cohort.path), current.cohort)) {
      stop("The cohort, reference, GVCFs, or filter settings changed since the ",
           "existing results in ", output.directory, ".\nRerun with overwrite = TRUE ",
           "to regenerate all loci, or use a new output directory.")
    }
  }

  # Creates every owned subdirectory recursively on each run
  final.dirs = c("filtered-all", "filtered-snps", "filtered-indels")
  work.dirs = c("unfiltered-all", "unfiltered-snps", "unfiltered-indels")
  workspace.root = file.path(output.directory, "genomicsdb-workspace")
  completion.root = file.path(output.directory, "completion")
  for (d in c(final.dirs, work.dirs, "completion")) {
    .ensureDirectory(file.path(output.directory, d))
  }
  .ensureDirectory(workspace.root)
  saveRDS(current.cohort, cohort.path)

  # Requested final products per locus, derived from the save flags
  requested = character(0)
  if (save.combined) requested = c(requested, "filtered-all")
  if (save.SNPs) requested = c(requested, "filtered-snps")
  if (save.indels) requested = c(requested, "filtered-indels")

  # One validated sample-name map is reused for every locus import, instead of a
  # long repeated -V argument string.
  sample.map = file.path(output.directory, "cohort-sample-map.txt")
  write.table(data.frame(vcf.sample.ids, normalizePath(gvcf.files)),
              sample.map, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Filter strings, unchanged thresholds
  SNP.QD.string <- paste0(" -filter \"QD<", format(custom.SNP.QD, nsmall = 1), "\" --filter-name \"QD\"")
  SNP.QUAL.string <- paste0(" -filter \"QUAL<", format(custom.SNP.QUAL, nsmall = 1), "\" --filter-name \"QUAL\"")
  SNP.SOR.string <- paste0(" -filter \"SOR>", format(custom.SNP.SOR, nsmall = 1), "\" --filter-name \"SOR\"")
  SNP.FS.string <- paste0(" -filter \"FS>", format(custom.SNP.FS, nsmall = 1), "\" --filter-name \"FS\"")
  SNP.MQ.string <- paste0(" -filter \"MQ<", format(custom.SNP.MQ, nsmall = 1), "\" --filter-name \"MQ\"")
  SNP.MQRankSum.string <- paste0(" -filter \"MQRankSum<", format(custom.SNP.MQRankSum, nsmall = 1), "\" --filter-name \"MQRankSum\"")
  SNP.ReadPosRankSum.string <- paste0(" -filter \"ReadPosRankSum<", format(custom.SNP.ReadPosRankSum, nsmall = 1), "\" --filter-name \"ReadPosRankSum\"")
  IN.QD.string <- paste0(" -filter \"QD<", format(custom.INDEL.QD, nsmall = 1), "\" --filter-name \"QD\"")
  IN.QUAL.string <- paste0(" -filter \"QUAL<", format(custom.INDEL.QUAL, nsmall = 1), "\" --filter-name \"QUAL\"")
  IN.FS.string <- paste0(" -filter \"FS>", format(custom.INDEL.FS, nsmall = 1), "\" --filter-name \"FS\"")
  IN.ReadPosRankSum.string <- paste0(" -filter \"ReadPosRankSum<", format(custom.INDEL.ReadPosRankSum, nsmall = 1), "\" --filter-name \"ReadPosRankSum\"")

  # Selects loci that still need work
  pending = loci.names[!vapply(loci.names, function(locus) {
    outputs = file.path(output.directory, requested, paste0(locus, ".vcf"))
    .stageComplete(file.path(completion.root, locus), "jointGenotyping", outputs)
  }, logical(1))]

  if (length(pending) == 0) {
    .cleanupJointOutputs(output.directory, save.unfiltered, save.SNPs, save.indels, save.combined)
    return(invisible(loci.names))
  }

  resources = .validateResources(threads, memory, length(pending))
  gatk = .toolCommand("gatk", gatk4.path)
  gatk.command = .gatkCommand(gatk, temp.directory, resources$heap.mb)

  results = parallel::mclapply(seq_along(pending), function(i) {
    locus = pending[i]
    log = file.path("logs/sample_logs", paste0("FAILURE_", locus, "_jointGenotyping.txt"))
    tryCatch({
      unfiltered.all = file.path(output.directory, "unfiltered-all", paste0(locus, ".vcf"))
      unfiltered.snps = file.path(output.directory, "unfiltered-snps", paste0(locus, ".vcf"))
      unfiltered.indels = file.path(output.directory, "unfiltered-indels", paste0(locus, ".vcf"))
      workspace = file.path(workspace.root, locus)

      # GenomicsDBImport requires an empty workspace path, so any leftover
      # workspace is removed before the import.
      if (dir.exists(workspace)) { unlink(workspace, recursive = TRUE) }

      .runCommand(paste0(gatk.command, " GenomicsDBImport",
                    " --sample-name-map ", shQuote(sample.map),
                    " --genomicsdb-workspace-path ", shQuote(workspace),
                    " --intervals ", shQuote(locus),
                    " --batch-size ", batch.size,
                    " --reader-threads 1"),
                  quiet, "GenomicsDBImport", stderr.log = log)

      .runCommand(paste0(gatk.command, " GenotypeGVCFs -R ", shQuote(reference.path),
                    " -V ", shQuote(paste0("gendb://", workspace)),
                    " --use-new-qual-calculator true",
                    " -O ", shQuote(unfiltered.all)),
                  quiet, "GenotypeGVCFs", stderr.log = log)

      .runCommand(paste0(gatk.command, " SelectVariants -V ", shQuote(unfiltered.all),
                    " -O ", shQuote(unfiltered.snps), " --select-type SNP"),
                  quiet, "SNP selection", stderr.log = log)
      .runCommand(paste0(gatk.command, " SelectVariants -V ", shQuote(unfiltered.all),
                    " -O ", shQuote(unfiltered.indels), " --select-type INDEL"),
                  quiet, "indel selection", stderr.log = log)

      snps.filter = sub("\\.vcf$", "_filter.vcf", unfiltered.snps)
      indels.filter = sub("\\.vcf$", "_filter.vcf", unfiltered.indels)
      all.filter = sub("\\.vcf$", "_filter.vcf", unfiltered.all)

      .runCommand(paste0(gatk.command, " VariantFiltration -R ", shQuote(reference.path),
                    " -V ", shQuote(unfiltered.snps), " -O ", shQuote(snps.filter),
                    SNP.QD.string, SNP.QUAL.string, SNP.SOR.string, SNP.FS.string,
                    SNP.MQ.string, SNP.MQRankSum.string, SNP.ReadPosRankSum.string),
                  quiet, "SNP filtering", stderr.log = log)
      .runCommand(paste0(gatk.command, " VariantFiltration -R ", shQuote(reference.path),
                    " -V ", shQuote(unfiltered.indels), " -O ", shQuote(indels.filter),
                    IN.QD.string, IN.QUAL.string, IN.FS.string, IN.ReadPosRankSum.string),
                  quiet, "indel filtering", stderr.log = log)

      .runCommand(paste0(gatk.command, " SortVcf -I ", shQuote(snps.filter),
                    " -I ", shQuote(indels.filter), " -O ", shQuote(all.filter)),
                  quiet, "VCF merge", stderr.log = log)

      # Publishes the requested final products
      if (save.combined) {
        .runCommand(paste0(gatk.command, " SelectVariants -V ", shQuote(all.filter),
                      " -O ", shQuote(file.path(output.directory, "filtered-all", paste0(locus, ".vcf"))),
                      " --exclude-filtered TRUE"),
                    quiet, "combined passing selection", stderr.log = log)
      }
      if (save.SNPs) {
        .runCommand(paste0(gatk.command, " SelectVariants -V ", shQuote(snps.filter),
                      " -O ", shQuote(file.path(output.directory, "filtered-snps", paste0(locus, ".vcf"))),
                      " --exclude-filtered TRUE"),
                    quiet, "SNP passing selection", stderr.log = log)
      }
      if (save.indels) {
        .runCommand(paste0(gatk.command, " SelectVariants -V ", shQuote(indels.filter),
                      " -O ", shQuote(file.path(output.directory, "filtered-indels", paste0(locus, ".vcf"))),
                      " --exclude-filtered TRUE"),
                    quiet, "indel passing selection", stderr.log = log)
      }

      outputs = file.path(output.directory, requested, paste0(locus, ".vcf"))
      if (!all(file.exists(outputs) & file.info(outputs)$size > 0)) {
        stop("Requested joint-genotyping products are incomplete for locus ", locus)
      }
      .ensureDirectory(file.path(completion.root, locus))
      .markStageComplete(file.path(completion.root, locus), "jointGenotyping",
                         paste0("cohort=", length(sample.names)))
      # Removes only this completed locus's owned workspace
      if (dir.exists(workspace)) { unlink(workspace, recursive = TRUE) }
      list(success = TRUE)
    }, error = function(e) {
      cat("\n", conditionMessage(e), "\n", file = log, append = TRUE)
      list(success = FALSE, message = conditionMessage(e))
    })
  }, mc.cores = resources$workers)

  .collectWorkers(results, pending, "Joint genotyping")

  # Optional cleanup runs only after every requested product succeeded
  .cleanupJointOutputs(output.directory, save.unfiltered, save.SNPs, save.indels, save.combined)
  invisible(loci.names)
}#end function


# Removes optional VCF subdirectories after successful publication of every
# requested final product.
.cleanupJointOutputs = function(output.directory, save.unfiltered, save.SNPs,
                                save.indels, save.combined) {
  if (save.SNPs == FALSE) {
    unlink(file.path(output.directory, c("unfiltered-snps", "filtered-snps")), recursive = TRUE)
  }
  if (save.indels == FALSE) {
    unlink(file.path(output.directory, c("unfiltered-indels", "filtered-indels")), recursive = TRUE)
  }
  if (save.combined == FALSE) {
    unlink(file.path(output.directory, c("unfiltered-all", "filtered-all")), recursive = TRUE)
  }
  if (save.unfiltered == FALSE) {
    unlink(file.path(output.directory, c("unfiltered-all", "unfiltered-snps", "unfiltered-indels")),
           recursive = TRUE)
  }
  invisible(NULL)
}


# Reads the header lines of a plain or gzipped VCF up to and including the
# #CHROM line. Avoids loading records into memory.
.vcfHeaderLines = function(path) {
  con = if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  header = character(0)
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) { header = c(header, line); break }
    if (startsWith(line, "##")) { header = c(header, line) }
  }
  header
}

# Sample IDs are the header columns after the fixed FORMAT column.
.vcfSamples = function(header.lines) {
  chrom = header.lines[startsWith(header.lines, "#CHROM")]
  if (length(chrom) != 1) { return(character(0)) }
  fields = strsplit(chrom, "\t")[[1]]
  if (length(fields) <= 9) { return(character(0)) }
  fields[10:length(fields)]
}

# Contig names declared in the header.
.vcfContigs = function(header.lines) {
  contig.lines = header.lines[startsWith(header.lines, "##contig=")]
  if (length(contig.lines) == 0) { return(character(0)) }
  sub(".*ID=([^,>]+).*", "\\1", contig.lines)
}

# END SCRIPT
