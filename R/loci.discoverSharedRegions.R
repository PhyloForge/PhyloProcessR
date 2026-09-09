#' @title discoverSharedRegions
#'
#' @description Discovers novel genomic loci by identifying regions of a reference genome
#' that are consistently covered by reads that did not map to any known sequence-capture locus.
#' For each sample, reads are first mapped to a consensus reference built from existing
#' capture alignments; unmapped reads are then mapped to the reference genome. Genomic
#' regions covered in at least \code{min.samples} samples at a depth of at least
#' \code{min.coverage} reads are retained as candidate novel loci. The genome sequence at
#' each shared region is extracted and written as a FASTA target file for use with downstream
#' annotation and alignment functions. Per-sample genome-mapped BAM files are retained for
#' \code{assembleSharedRegions}.
#'
#' @param alignment.directory path to the directory containing existing sequence-capture
#' alignment files (e.g. \code{data-analysis/alignments/untrimmed_all-markers}).
#'
#' @param alignment.format format of the alignment files. Accepted: "phylip" or "fasta".
#' Default "phylip".
#'
#' @param read.directory path to the directory containing per-sample read subdirectories.
#' Point this directly at the folder whose immediate children are one directory per sample
#' (e.g. \code{"processed-reads/decontaminated-reads"}). R1/R2 files are expected directly
#' inside each sample subdirectory.
#'
#' @param genome.file full path to the reference genome FASTA file.
#'
#' @param output.directory path to write output files (BAMs, BEDs, target FASTA).
#'
#' @param min.samples minimum number of samples that must cover a region for it to be
#' retained as a candidate novel locus. Default 4.
#'
#' @param min.coverage minimum read depth required at a site for it to be considered
#' covered in a sample. Default 5.
#'
#' @param min.region.length minimum length in bp for a candidate novel region to be retained.
#' Default 200.
#'
#' @param max.merge.distance maximum distance in bp between adjacent covered intervals after
#' the shared-coverage threshold is applied. Gaps can be included in the final region. Default 500.
#'
#' @param min.mapping.quality minimum MAPQ score for a read to be counted. Filters
#' multi-mapping reads in repetitive regions. Default 20.
#'
#' @param threads number of CPU threads for parallel sample processing. Default 1.
#'
#' @param memory total RAM in GB available. Default 8.
#'
#' @param overwrite logical. If TRUE, rebuilds references and per-sample mapping outputs.
#' Use TRUE after changing inputs or mapping/coverage settings. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses stdout/stderr from external tools. Default FALSE.
#'
#' @param hisat2.path path to directory containing HISAT2 executables. Default NULL (system PATH).
#'
#' @param samtools.path path to directory containing samtools executable. Default NULL (system PATH).
#'
#' @param bedtools.path path to directory containing bedtools executable. Default NULL (system PATH).
#'
#' @return Writes to \code{output.directory}:
#' \itemize{
#'   \item \code{sample-bams/} -- per-sample genome-mapped BAM files
#'   \item \code{novel_regions.bed} -- BED file of shared novel regions
#'   \item \code{novel_targets.fa} -- genome sequences at shared regions (use as target file
#'         for \code{annotateTargets})
#' }
#' No value is returned to R.
#'
#' @export

discoverSharedRegions = function(alignment.directory = NULL,
                                 alignment.format = "phylip",
                                 read.directory = NULL,
                                 genome.file = NULL,
                                 output.directory = NULL,
                                 min.samples = 4,
                                 min.coverage = 5,
                                 min.region.length = 200,
                                 max.merge.distance = 500,
                                 min.mapping.quality = 20,
                                 threads = 1,
                                 memory = 8,
                                 overwrite = FALSE,
                                 quiet = FALSE,
                                 hisat2.path = NULL,
                                 samtools.path = NULL,
                                 bedtools.path = NULL) {

  # alignment.directory = "data-analysis/alignments/untrimmed_all-markers"
  # alignment.format = "phylip"
  # read.directory = "processed-reads/decontaminated-reads"
  # genome.file = "/PATH/TO/genome.fa"
  # output.directory = "data-analysis/novel-loci-discovery"
  # min.samples = 4
  # min.coverage = 5
  # min.region.length = 200
  # max.merge.distance = 500
  # min.mapping.quality = 20
  # threads = 8
  # memory = 40
  # overwrite = TRUE
  # quiet = TRUE
  # hisat2.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"
  # samtools.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"
  # bedtools.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"

  #Input checks
  if (is.null(alignment.directory)) stop("alignment.directory not provided.")
  if (is.null(read.directory)) stop("read.directory not provided.")
  if (is.null(genome.file)) stop("genome.file not provided.")
  if (is.null(output.directory)) stop("output.directory not provided.")
  if (!file.exists(genome.file)) {
    stop("genome.file not found: ", genome.file)
  }

  hisat2.command = .toolCommand("hisat2", hisat2.path)
  build.command = .toolCommand("hisat2-build", hisat2.path)
  samtools.command = .toolCommand("samtools", samtools.path)
  bedtools.command = .toolCommand("bedtools", bedtools.path)

  # Create output directories; never wipe on resume -- overwrite controls per-step redo
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(output.directory, "/sample-bams"), showWarnings = FALSE)
  dir.create(paste0(output.directory, "/covered-beds"), showWarnings = FALSE)

  #Gather alignments
  if (!alignment.format %in% c("phylip", "fasta")) {
    stop("alignment.format must be phylip or fasta.")
  }
  if (alignment.format == "phylip") {
    alignment.pattern = "\\.(phy|phylip)$"
  } else {
    alignment.pattern = "\\.(fa|fas|fasta|fna)$"
  }
  align.files = list.files(alignment.directory, pattern = alignment.pattern,
                           ignore.case = TRUE)
  if (length(align.files) == 0) stop("No alignment files found in alignment.directory.")

  known.ref    = paste0(output.directory, "/known_loci_consensus.fa")
  known.index  = paste0(output.directory, "/known_loci_index")
  genome.index = paste0(output.directory, "/genome_index")

  ##################################################################################################
  ## Step 1: Build known-loci consensus reference
  ##################################################################################################
  rebuild.known = overwrite || !file.exists(known.ref) || file.size(known.ref) == 0
  if (rebuild.known) {
    print(paste0("Building known-loci consensus from ", length(align.files), " alignments..."))

    all.consensus = Biostrings::DNAStringSet()
    for (i in 1:length(align.files)) {
      locus.name = tools::file_path_sans_ext(align.files[i])
      if (alignment.format == "phylip") {
        align = Biostrings::DNAStringSet(Biostrings::readDNAMultipleAlignment(
          file = paste0(alignment.directory, "/", align.files[i]), format = "phylip"))
      } else {
        align = Biostrings::readDNAStringSet(paste0(alignment.directory, "/", align.files[i]))
      }
      con = makeConsensus(alignment = align, method = "majority",
                          warn.non.IUPAC = FALSE, remove.gaps = TRUE, type = "DNA")
      names(con) = locus.name
      all.consensus = append(all.consensus, con)
      rm(align, con)
    }

    Biostrings::writeXStringSet(all.consensus, filepath = known.ref)
    rm(all.consensus)
    gc()
  } else {
    print("Known-loci consensus already exists -- skipping Step 1.")
  }

  ##################################################################################################
  ## Step 2: Index references
  ##################################################################################################
  index.prefixes = c(known.index, genome.index)
  reference.files = c(known.ref, genome.file)
  for (i in seq_along(index.prefixes)) {
    small.index = paste0(index.prefixes[i], ".", 1:8, ".ht2")
    large.index = paste0(index.prefixes[i], ".", 1:8, ".ht2l")
    small.complete = all(file.exists(small.index)) && all(file.size(small.index) > 0)
    large.complete = all(file.exists(large.index)) && all(file.size(large.index) > 0)
    if (overwrite || (i == 1 && rebuild.known) || !(small.complete || large.complete)) {
      unlink(c(small.index, large.index))
      .runCommand(paste(build.command, shQuote(reference.files[i]),
                         shQuote(index.prefixes[i])),
                  quiet = quiet, task = "HISAT2 indexing")
    }
  }

  ##################################################################################################
  ## Step 3: Per-sample -- map to known loci, extract unmapped, map to genome
  ##################################################################################################
  sample.names = list.dirs(read.directory, recursive = FALSE, full.names = FALSE)
  if (length(sample.names) == 0) {
    stop("No sample directories found in read.directory.")
  }
  print(paste0("Processing ", length(sample.names), " samples..."))

  sample.results = parallel::mclapply(sample.names, function(samp) {
    bed.out = file.path(output.directory, "covered-beds", paste0(samp, "_covered.bed"))
    bam.out = file.path(output.directory, "sample-bams", paste0(samp, ".bam"))
    done.file = paste0(bed.out, ".complete")
    tryCatch({
      if (!overwrite && file.exists(done.file) && file.exists(bed.out) &&
          file.exists(bam.out) && file.exists(paste0(bam.out, ".bai"))) {
        return(TRUE)
      }
      unlink(done.file)
      read.dir = file.path(read.directory, samp)
      read.files = .listFastqFiles(read.dir, recursive = FALSE)
      lane.prefixes = .stripReadSuffix(read.files)
      if (length(lane.prefixes) == 0) stop("No paired FASTQ files found.")
      pairs = lapply(lane.prefixes, function(prefix) {
        lane.files = .matchPrefix(read.files, read.files, prefix)
        .orderReadFiles(lane.files)
      })
      r1 = vapply(pairs, function(pair) pair[1], character(1))
      r2 = vapply(pairs, function(pair) pair[2], character(1))

      tmp = file.path(output.directory, paste0("tmp_", samp))
      dir.create(tmp, showWarnings = FALSE)
      known.sam = file.path(tmp, "known.sam")
      unmapped.bam = file.path(tmp, "unmapped_sorted.bam")
      unmapped.r1 = file.path(tmp, "unmapped_R1.fastq")
      unmapped.r2 = file.path(tmp, "unmapped_R2.fastq")
      genome.sam = file.path(tmp, "genome.sam")

      # Keep pairs for which neither mate maps to the known loci.
      .runCommand(paste0(hisat2.command, " -x ", shQuote(known.index),
                          " -1 ", shQuote(paste(r1, collapse = ",")),
                          " -2 ", shQuote(paste(r2, collapse = ",")),
                          " --threads 1 --no-spliced-alignment -S ", shQuote(known.sam)),
                  quiet = quiet, task = "known-locus mapping")
      .runPipeline(paste0(samtools.command, " view -b -f 12 -F 2304 ", shQuote(known.sam),
                           " | ", samtools.command, " sort -n -o ", shQuote(unmapped.bam)),
                   quiet = quiet, task = "unmapped read extraction")
      .runCommand(paste0(samtools.command, " fastq -1 ", shQuote(unmapped.r1),
                          " -2 ", shQuote(unmapped.r2),
                          " -0 /dev/null -s /dev/null -n ", shQuote(unmapped.bam)),
                  quiet = quiet, task = "unmapped FASTQ conversion")

      # HISAT2 accepts empty FASTQs and writes a header-only SAM when no pairs remain.
      .runCommand(paste0(hisat2.command, " -x ", shQuote(genome.index),
                          " -1 ", shQuote(unmapped.r1), " -2 ", shQuote(unmapped.r2),
                          " --threads 1 --no-spliced-alignment -S ", shQuote(genome.sam)),
                  quiet = quiet, task = "genome mapping")
      .runPipeline(paste0(samtools.command, " view -b -F 2308 -q ", min.mapping.quality,
                           " ", shQuote(genome.sam), " | ", samtools.command,
                           " sort -o ", shQuote(bam.out)),
                   quiet = quiet, task = "genome BAM sorting")
      .runCommand(paste(samtools.command, "index", shQuote(bam.out)),
                  quiet = quiet, task = "genome BAM indexing")

      # Keep actual covered bases. Gap merging occurs after the sample threshold.
      bed.temp = paste0(bed.out, ".tmp")
      .runPipeline(paste0(bedtools.command, " genomecov -ibam ", shQuote(bam.out),
                           " -bg | awk '$4 >= ", min.coverage, "' | ",
                           bedtools.command, " merge -i stdin > ", shQuote(bed.temp)),
                   quiet = quiet, keep.stdout = TRUE, task = "coverage intervals")
      if (!file.rename(bed.temp, bed.out)) stop("Cannot save coverage BED.")
      file.create(done.file)
      unlink(tmp, recursive = TRUE)
      print(paste0("Finished mapping ", samp))
      TRUE
    }, error = function(e) {
      message("Error processing ", samp, ": ", conditionMessage(e))
      FALSE
    })
  }, mc.cores = threads)
  if (!all(vapply(sample.results, isTRUE, logical(1)))) {
    stop("Discovery failed for one or more samples. Correct the errors and rerun.")
  }

  ##################################################################################################
  ## Step 4: Find genomic regions shared across >= min.samples samples
  ##################################################################################################
  shared.bed = paste0(output.directory, "/novel_regions.bed")
  novel.fa   = paste0(output.directory, "/novel_targets.fa")

  writeLines(sample.names, file.path(output.directory, "samples.txt"))

  # Use only the samples in this run. Recompute shared outputs from their BEDs.
  bed.files = file.path(output.directory, "covered-beds",
                         paste0(sample.names, "_covered.bed"))
  .sharedCoveredRegions(bed.files, shared.bed, min.samples, min.region.length,
                         max.merge.distance, bedtools.command, quiet)
  if (file.size(shared.bed) == 0) {
    writeLines(character(), novel.fa)
    print("No shared novel regions passed the coverage and length thresholds.")
    return(invisible(NULL))
  }

  ##################################################################################################
  ## Step 5: Extract genome sequences at shared regions
  ##################################################################################################
  raw.fa = paste0(output.directory, "/novel_targets_raw.fa")
  .runCommand(paste0(bedtools.command, " getfasta -fi ", shQuote(genome.file),
                      " -bed ", shQuote(shared.bed), " -fo ", shQuote(raw.fa)),
              quiet = quiet, task = "novel target extraction")
  raw.seqs = Biostrings::readDNAStringSet(raw.fa)
  names(raw.seqs) = gsub("[^A-Za-z0-9_.]", "_", names(raw.seqs))
  if (anyDuplicated(names(raw.seqs))) {
    stop("Genome scaffold names produce duplicate novel target names after sanitization.")
  }
  Biostrings::writeXStringSet(raw.seqs, filepath = novel.fa)
  unlink(raw.fa)

  print(paste0("Novel target sequences written to: ", novel.fa))
  print(paste0("Per-sample BAM files written to:   ", output.directory, "/sample-bams/"))
  print(paste0("Shared regions BED:                ", shared.bed))

}#end function

#END SCRIPT
