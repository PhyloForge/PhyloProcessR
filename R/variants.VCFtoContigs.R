#' Apply filtered variants and depth rules to sample assemblies
#'
#' Produces an alternate-reference sequence from passing variants, optionally
#' encoding supported heterozygous SNP genotypes with IUPAC codes. Reference
#' bases remain outside applied variants. This is not phased haplotype output.
#'
#' @param genotype.directory Per-sample genotype directory.
#' @param mapping.directory Per-sample mapping/reference directory.
#' @param output.directory FASTA output directory. When both sequence modes are
#'   selected, this is the parent of `4_consensus-contigs` and `5_iupac-contigs`.
#' @param vcf.file One of `SNP`, `Indel`, or `Both`.
#' @param consensus.sequences Produce ordinary alternate-reference output.
#' @param ambiguity.codes Produce IUPAC output.
#' @param threads Concurrent samples.
#' @param memory Total JVM heap budget in GB.
#' @param temp.directory Temporary directory.
#' @param gatk4.path GATK path/directory or NULL.
#' @param overwrite Recompute output.
#' @param quiet Suppress tool output while retaining logs.
#' @param sample.names Optional retained sample set.
#' @param depth.files Optional named depth-table paths shared between output modes.
#' @param depth.filter.mode `none`, `site`, `mean`, or `both` (default `site`).
#' @param min.site.depth Sites below this depth become N; default 1 retains 1x.
#' @param min.mean.depth Contigs below this full-span mean are removed.
#' @param max.n.proportion Optional final N-fraction cutoff.
#' @param use.base.recalibration Select recalibrated BAM when calculating depth.
#' @param samtools.path samtools path/directory or NULL.
#' @param ploidy Used only for a targeted non-diploid IUPAC warning.
#' @return Invisibly returns sample names.
#' @export
VCFtoContigs = function(genotype.directory = NULL,
                        mapping.directory = NULL,
                        output.directory = "sample-contigs",
                        vcf.file = "SNP",
                        consensus.sequences = FALSE,
                        ambiguity.codes = TRUE,
                        threads = 1,
                        memory = 1,
                        temp.directory = NULL,
                        gatk4.path = NULL,
                        overwrite = FALSE,
                        quiet = TRUE,
                        sample.names = NULL,
                        depth.files = NULL,
                        depth.filter.mode = "site",
                        min.site.depth = 1,
                        min.mean.depth = 1,
                        max.n.proportion = NULL,
                        use.base.recalibration = FALSE,
                        samtools.path = NULL,
                        ploidy = 2) {

  #################################################
  ### Part A: check the arguments and inputs
  #################################################
  choice = match.arg(tolower(vcf.file), c("snp", "indel", "both"))
  vcf.file = c(snp = "SNP", indel = "Indel", both = "Both")[[choice]]

  if (!isTRUE(consensus.sequences) && !isTRUE(ambiguity.codes)) {
    stop("At least one of consensus.sequences and ambiguity.codes must be TRUE.")
  }

  settings = .validateDepthSettings(depth.filter.mode, min.site.depth,
                                    min.mean.depth, max.n.proportion)
  if (vcf.file != "SNP" && settings$mode %in% c("site", "both")) {
    stop("Site depth masking currently requires vcf.file = 'SNP'; ",
         "use mean or none for indel output.")
  }
  if (ambiguity.codes == TRUE && ploidy != 2) {
    warning("GATK IUPAC output does not fully encode non-diploid genotypes.")
  }

  if (is.null(genotype.directory) || !dir.exists(genotype.directory)) {
    stop("Genotype directory not found.")
  }
  if (is.null(mapping.directory) || !dir.exists(mapping.directory)) {
    stop("Mapping directory not found.")
  }

  if (is.null(sample.names)) {
    sample.names = list.dirs(genotype.directory, recursive = FALSE, full.names = FALSE)
  }
  if (length(sample.names) == 0) {
    stop("No samples are available for FASTA conversion.")
  }

  #Selects the VCF for the requested variant type and the per-sample reference
  vcf.name = switch(vcf.file,
                    SNP = "gatk4-final-snps.vcf",
                    Indel = "gatk4-final-indels.vcf",
                    Both = "gatk4-final-genotypes.vcf")
  vcf.paths = file.path(genotype.directory, sample.names, vcf.name)
  ref.paths = file.path(mapping.directory, sample.names, "index", "reference.fa")
  if (any(!file.exists(vcf.paths))) {
    stop("Missing selected VCF for sample(s): ",
         paste(sample.names[!file.exists(vcf.paths)], collapse = ", "))
  }
  if (any(!file.exists(ref.paths))) {
    stop("Missing reference for sample(s): ",
         paste(sample.names[!file.exists(ref.paths)], collapse = ", "))
  }

  #Write both formats to separate directories. Reuse the same depth tables so
  #the two outputs apply the same depth rules to each sample.
  if (isTRUE(consensus.sequences) && isTRUE(ambiguity.codes)) {
    need.depth = settings$mode != "none" || !is.null(max.n.proportion)
    if (need.depth && is.null(depth.files)) {
      depth.files = calculateSampleDepth(mapping.directory,
                                         file.path(output.directory, "depth"),
                                         sample.names, use.base.recalibration,
                                         samtools.path, overwrite, quiet)
    }
    for (mode in c("consensus", "iupac")) {
      VCFtoContigs(
        genotype.directory = genotype.directory,
        mapping.directory = mapping.directory,
        output.directory = file.path(output.directory,
                                     if (mode == "consensus") "4_consensus-contigs" else "5_iupac-contigs"),
        vcf.file = vcf.file,
        consensus.sequences = mode == "consensus",
        ambiguity.codes = mode == "iupac",
        threads = threads,
        memory = memory,
        temp.directory = temp.directory,
        gatk4.path = gatk4.path,
        overwrite = overwrite,
        quiet = quiet,
        sample.names = sample.names,
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
    return(invisible(sample.names))
  }

  #################################################
  ### Part B: prepare depth tables and output paths
  #################################################
  #Depth is needed for site/mean masking or for the optional final N cutoff
  need.depth = settings$mode != "none" || !is.null(max.n.proportion)
  if (need.depth == TRUE && is.null(depth.files)) {
    depth.files = calculateSampleDepth(mapping.directory,
                                       file.path(dirname(output.directory), "depth"),
                                       sample.names, use.base.recalibration,
                                       samtools.path, overwrite, quiet)
  }

  resources = .validateResources(threads, memory, length(sample.names))
  gatk = .toolCommand("gatk", gatk4.path)

  if (is.null(temp.directory)) { temp.directory = tempdir() }
  .ensureDirectory(temp.directory)

  if (overwrite == TRUE && dir.exists(output.directory)) {
    unlink(output.directory, recursive = TRUE)
  }
  .ensureDirectory(output.directory)
  .ensureDirectory("logs/sample_logs")

  command = .gatkCommand(gatk, temp.directory, resources$heap.mb)
  mode = if (ambiguity.codes == TRUE) "iupac" else "consensus"

  #################################################
  ### Part C: convert each sample in parallel
  #################################################
  results = parallel::mclapply(seq_along(sample.names), function(i) {
    tryCatch({
      sample = sample.names[i]
      final.file = file.path(output.directory, paste0(sample, ".fa"))
      marker = .stageMarker(output.directory, paste0(sample, "-VCFtoContigs-", mode))
      report = file.path(output.directory, paste0(sample, ".", mode, ".depth-filter.tsv"))

      #The completion marker records the inputs and settings, so a rerun with the
      #same settings skips the sample and a changed setting recomputes it.
      details = c(paste0("vcf=", normalizePath(vcf.paths[i])),
                  paste0("mode=", settings$mode),
                  paste0("min.site.depth=", settings$min.site.depth),
                  paste0("min.mean.depth=", settings$min.mean.depth),
                  paste0("max.n.proportion=",
                         if (is.null(settings$max.n.proportion)) "NULL" else settings$max.n.proportion))
      if (overwrite == FALSE && file.exists(marker) && file.exists(final.file) &&
          identical(readLines(marker, warn = FALSE), c("complete=true", details))) {
        return(list(success = TRUE))
      }
      file.remove(marker)

      #Writes to temporary files first, then renames, so an interrupted run
      #cannot leave a partial FASTA at the final path.
      temp.raw = tempfile(paste0(sample, "-raw-"), tmpdir = output.directory, fileext = ".fa")
      temp.done = tempfile(paste0(sample, "-processed-"), tmpdir = output.directory, fileext = ".fa")
      log = file.path("logs/sample_logs", sample, paste0("VCFtoContigs-", mode, ".stderr.log"))
      .ensureDirectory(dirname(log))

      iupac = if (ambiguity.codes == TRUE) paste("--use-iupac-sample", shQuote(sample)) else ""
      .runCommand(paste(command, "FastaAlternateReferenceMaker -R", shQuote(ref.paths[i]),
                        "-V", shQuote(vcf.paths[i]), "-O", shQuote(temp.raw), iupac),
                  quiet, "FASTA conversion", stderr.log = log)
      if (!file.exists(temp.raw)) { stop("GATK did not create its temporary FASTA") }

      #GATK adds an index and full-contig interval to each reference name.
      sequences = Biostrings::readDNAStringSet(temp.raw)
      names(sequences) = .referenceContigNames(sequences)

      if (need.depth == TRUE) {
        depth.file = depth.files[[sample]]
        if (is.null(depth.file) || !file.exists(depth.file)) {
          stop("Depth result missing for ", sample)
        }
        sequences = .filterDepthSequences(sequences, depth.file, settings, report)
      }

      Biostrings::writeXStringSet(sequences, temp.done, format = "fasta", width = 1000000)
      if (!file.rename(temp.done, final.file)) {
        stop("Could not publish final FASTA for ", sample)
      }
      file.remove(c(temp.raw, paste0(temp.raw, ".fai"), sub("\\.fa$", ".dict", temp.raw)))

      writeLines(c("complete=true", details), marker)
      list(success = TRUE)
    }, error = function(e) list(success = FALSE, message = conditionMessage(e)))
  }, mc.cores = resources$workers)

  .collectWorkers(results, sample.names, "FASTA conversion")
  invisible(sample.names)
}#end function
