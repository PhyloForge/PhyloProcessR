#' @title assembleSharedRegions
#'
#' @description Assembles per-sample contigs for each novel genomic region identified by
#' \code{discoverSharedRegions}. For each sample all reads overlapping any novel region are
#' extracted in a single BAM pass, assembled together in one SPAdes run, and the resulting
#' contigs are assigned to individual regions by BLASTing against the \code{novel_targets.fa}
#' produced by \code{discoverSharedRegions}. Only contigs assigned to a region that had
#' \code{min.reads.assemble} or more reads are retained. Contigs are named
#' \code{region_contig_N} (e.g. \code{chr3_450000_450800_contig_1}) and written as one
#' FASTA per sample to \code{output.directory}, compatible with the downstream
#' \code{filterHeterozygosity} -> \code{collectNovelContigs} -> \code{alignTargets} pipeline
#' (workflow X4).
#'
#' @param discover.directory path to the output directory from \code{discoverSharedRegions}.
#' Must contain \code{sample-bams/}, \code{novel_regions.bed}, and \code{novel_targets.fa}.
#'
#' @param output.directory path to write per-sample contig FASTA files
#' (e.g. \code{data-analysis/contigs/9_genome-contigs}).
#'
#' @param min.reads.assemble minimum number of reads mapping to a region for contigs
#' assembled to that region to be retained. Default 5.
#'
#' @param kmer.values integer vector of k-mer sizes passed to SPAdes. Default
#' \code{c(33, 55, 77, 99, 127)}.
#'
#' @param threads number of CPU threads for parallel sample processing. Default 1.
#'
#' @param memory total RAM in GB available. Default 8.
#'
#' @param overwrite logical. If TRUE, existing per-sample FASTAs are re-generated. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses stdout/stderr from external tools. Default FALSE.
#'
#' @param spades.path path to directory containing the SPAdes executable. Default NULL
#' (system PATH).
#'
#' @param samtools.path path to directory containing samtools executable. Default NULL
#' (system PATH).
#'
#' @param bedtools.path path to directory containing the bedtools executable. Default NULL
#' (system PATH).
#'
#' @param blast.path path to directory containing blastn and makeblastdb executables.
#' Default NULL (system PATH).
#'
#' @return Writes one FASTA file per sample to \code{output.directory}. Intermediate
#' working directories are written under \code{discover.directory} and cleaned up after
#' each sample completes. No value is returned to R.
#'
#' @export

assembleSharedRegions = function(discover.directory = NULL,
                                 output.directory = NULL,
                                 min.reads.assemble = 5,
                                 kmer.values = c(33, 55, 77, 99, 127),
                                 threads = 1,
                                 memory = 8,
                                 overwrite = FALSE,
                                 quiet = FALSE,
                                 spades.path = NULL,
                                 samtools.path = NULL,
                                 bedtools.path = NULL,
                                 blast.path = NULL) {

  # discover.directory = "data-analysis/novel-loci-discovery"
  # output.directory = "data-analysis/contigs/9_genome-contigs"
  # min.reads.assemble = 5
  # kmer.values = c(33, 55, 77, 99, 127)
  # threads = 8
  # memory = 40
  # overwrite = TRUE
  # quiet = TRUE
  # spades.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"
  # samtools.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"
  # bedtools.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"
  # blast.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"

  # Input checks
  if (is.null(discover.directory)) { print("discover.directory not provided."); return(NULL) }
  if (is.null(output.directory))   { print("output.directory not provided.");   return(NULL) }

  region.bed  = paste0(discover.directory, "/novel_regions.bed")
  novel.fa    = paste0(discover.directory, "/novel_targets.fa")
  bam.dir     = paste0(discover.directory, "/sample-bams")

  if (!file.exists(region.bed)) {
    stop("novel_regions.bed not found in: ", discover.directory)
  }
  if (!file.exists(novel.fa)) {
    stop("novel_targets.fa not found in: ", discover.directory)
  }
  if (!dir.exists(bam.dir)) {
    stop("sample-bams/ not found in: ", discover.directory)
  }

  spades.command = .toolCommand("spades.py", spades.path)
  samtools.command = .toolCommand("samtools", samtools.path)
  bedtools.command = .toolCommand("bedtools", bedtools.path)
  blast.command = .toolCommand("blastn", blast.path)
  makeblastdb.command = .toolCommand("makeblastdb", blast.path)
  if (memory < 4) stop("Assembly requires at least 4 GB of memory.")

  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  if (overwrite) {
    old.fastas = list.files(output.directory, pattern = "\\.fa$", full.names = TRUE)
    unlink(old.fastas)
  }

  # Load region BED
  if (file.size(region.bed) == 0) {
    print("novel_regions.bed is empty -- no shared regions were found.")
    return(invisible(NULL))
  }
  regions = data.table::fread(region.bed, sep = "\t", header = FALSE)
  if (is.null(regions) || nrow(regions) == 0) {
    print("novel_regions.bed is empty -- no shared regions were found.")
    return(NULL)
  }
  if (ncol(regions) != 4) stop("novel_regions.bed must contain four columns.")
  data.table::setnames(regions, c("chrom", "start", "end", "sample_count"))
  # Sanitize region names the same way discoverSharedRegions sanitizes novel_targets.fa
  # sequence names -- must match exactly so the BLAST tName lookup works.
  # Scaffold names from some genome assemblies contain shell-unsafe characters
  # (e.g. "ScWFrIx_100724;HRSCAF=196187") that break file paths and MAFFT commands.
  region.names = gsub("[^A-Za-z0-9_.]", "_",
                      paste0(regions$chrom, "_", regions$start, "_", regions$end))
  print(paste0("Assembling contigs for ", nrow(regions), " regions..."))

  # Get per-sample BAM files
  bam.files = list.files(bam.dir, pattern = "\\.bam$", full.names = TRUE)
  bam.files = bam.files[!grepl("\\.bai$", bam.files)]
  sample.file = file.path(discover.directory, "samples.txt")
  if (file.exists(sample.file)) {
    current.samples = readLines(sample.file)
    bam.files = file.path(bam.dir, paste0(current.samples, ".bam"))
    if (!all(file.exists(bam.files))) stop("A current sample BAM is missing.")
  }
  sample.names = gsub("\\.bam$", "", basename(bam.files))

  if (length(bam.files) == 0) stop("No BAM files found in sample-bams/.")
  print(paste0("Assembling ", length(sample.names), " samples..."))

  ##################################################################################################
  ## Build a shared BLAST database from novel_targets.fa (once, before parallel loop)
  ##################################################################################################
  blast.db = paste0(discover.directory, "/novel_targets_blast_db")
  .runCommand(paste0(makeblastdb.command, " -in ", shQuote(novel.fa),
                " -dbtype nucl -out ", shQuote(blast.db)),
         quiet = quiet)

  kmer.str  = paste0(kmer.values, collapse = ",")
  workers = min(threads, length(sample.names), floor(memory / 4))
  mem.per = floor(memory / workers)
  blast.headers = c("qName", "tName", "pident", "length", "mismatch", "gapopen",
                    "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")

  ##################################################################################################
  ## Per-sample assembly (parallel)
  ##################################################################################################
  sample.results = parallel::mclapply(seq_along(sample.names), function(s) {
    samp = sample.names[s]
    out.fa = file.path(output.directory, paste0(samp, ".fa"))
    done.file = file.path(output.directory, paste0(".", samp, ".complete"))
    tryCatch({

      bam  = bam.files[s]

      if (!overwrite && file.exists(out.fa) && file.exists(done.file)) {
        print(paste0(samp, ": contig FASTA already exists -- skipping."))
        return(TRUE)
      }

      unlink(c(out.fa, done.file))
      samp.dir = tempfile(paste0("assembly_", samp, "_"), tmpdir = discover.directory)
      dir.create(samp.dir, showWarnings = FALSE)

      ##########################################################################
      # Step A: Extract all reads overlapping any novel region in one BAM pass
      ##########################################################################
      novel.bam = paste0(samp.dir, "/novel_reads.bam")
      ret0 = .runCommand(paste0(samtools.command, " view -b -L ", shQuote(region.bed),
                           " -o ", shQuote(novel.bam), " ", shQuote(bam)),
                    quiet = quiet)
      if (ret0 != 0 || !file.exists(novel.bam) || file.size(novel.bam) == 0) {
        unlink(samp.dir, recursive = TRUE)
        stop("Failed to extract novel-region reads for ", samp)
      }
      .runCommand(paste0(samtools.command, " index ", shQuote(novel.bam)),
             quiet = quiet)

      n.novel = as.integer(trimws(
        .runCommandOutput(paste0(samtools.command, " view -c ", shQuote(novel.bam)))))
      print(paste0(samp, ": ", n.novel, " reads mapped to novel regions."))

      if (is.na(n.novel) || n.novel == 0) {
        unlink(samp.dir, recursive = TRUE)
        print(paste0(samp, ": no reads in novel regions -- skipping."))
        writeLines(character(), out.fa)
        file.create(done.file)
        return(TRUE)
      }

      ##########################################################################
      # Step B: Count reads per region in one bedtools pass to determine which
      # regions have sufficient coverage (used to filter contigs after BLAST)
      ##########################################################################
      cov.file = file.path(samp.dir, "region_counts.bed")
      .runCommand(paste0(bedtools.command, " coverage -a ", shQuote(region.bed),
                          " -b ", shQuote(novel.bam), " -counts > ", shQuote(cov.file)),
                  quiet = quiet, keep.stdout = TRUE, task = "per-region read counts")
      region.counts = data.table::fread(cov.file, header = FALSE)
      if (nrow(region.counts) != nrow(regions) || ncol(region.counts) != 5) {
        stop("Unexpected region count table for ", samp)
      }
      reads.per.region = region.counts[[ncol(region.counts)]]
      active.regions   = region.names[reads.per.region >= min.reads.assemble]

      print(paste0(samp, ": ", length(active.regions), " of ", nrow(regions),
                   " regions have >= ", min.reads.assemble, " reads."))

      if (length(active.regions) == 0) {
        unlink(samp.dir, recursive = TRUE)
        print(paste0(samp, ": no regions with sufficient reads -- skipping."))
        writeLines(character(), out.fa)
        file.create(done.file)
        return(TRUE)
      }

      ##########################################################################
      # Step C: Convert all novel reads to FASTQ and run ONE SPAdes assembly
      ##########################################################################
      novel.fq   = paste0(samp.dir, "/novel_reads.fastq")
      spades.dir = paste0(samp.dir, "/spades_all")

      .runPipeline(paste0(samtools.command, " collate -u -O ", shQuote(novel.bam),
                           " | ", samtools.command, " fastq -N - > ", shQuote(novel.fq)),
                   quiet = quiet, keep.stdout = TRUE, task = "assembly FASTQ conversion")
      unlink(c(novel.bam, paste0(novel.bam, ".bai")))

      .runCommand(paste0(spades.command, " -s ", shQuote(novel.fq),
                    " -k ", kmer.str,
                    " -o ", shQuote(spades.dir),
                    " -m ", mem.per,
                    " --threads 1 --careful"),
             quiet = quiet)
      unlink(novel.fq)

      contig.fa = paste0(spades.dir, "/contigs.fasta")
      if (!file.exists(contig.fa) || file.size(contig.fa) == 0) {
        unlink(samp.dir, recursive = TRUE)
        print(paste0(samp, ": SPAdes produced no contigs."))
        writeLines(character(), out.fa)
        file.create(done.file)
        return(TRUE)
      }

      ##########################################################################
      # Step D: BLAST contigs against novel_targets.fa to assign region names.
      # tName values in the output are chr_start_end, matching region.names.
      ##########################################################################
      blast.out = paste0(samp.dir, "/contig_blast.txt")
      .runCommand(paste0(blast.command, " -task blastn -db ", shQuote(blast.db),
                    " -query ", shQuote(contig.fa),
                    " -out ", shQuote(blast.out),
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen",
                    " qstart qend sstart send evalue bitscore\"",
                    " -num_threads 1 -evalue 0.001"),
             quiet = quiet)

      ctgs = Biostrings::readDNAStringSet(contig.fa)
      unlink(spades.dir, recursive = TRUE)

      if (!file.exists(blast.out) || file.size(blast.out) == 0) {
        unlink(samp.dir, recursive = TRUE)
        print(paste0(samp, ": no BLAST hits for assembled contigs."))
        writeLines(character(), out.fa)
        file.create(done.file)
        return(TRUE)
      }

      blast.data = data.table::fread(blast.out, header = FALSE)
      unlink(blast.out)

      if (nrow(blast.data) == 0) {
        unlink(samp.dir, recursive = TRUE)
        print(paste0(samp, ": no BLAST hits for assembled contigs."))
        writeLines(character(), out.fa)
        file.create(done.file)
        return(TRUE)
      }
      data.table::setnames(blast.data, blast.headers)

      ##########################################################################
      # Step E: Assign each contig to its best-matching region, keep only
      # contigs assigned to regions that pass the min.reads.assemble threshold
      ##########################################################################
      # Best hit per contig = highest bitscore
      blast.best = blast.data[, .SD[which.max(bitscore)], by = qName]
      # Restrict to regions with sufficient read depth
      blast.best = blast.best[tName %in% active.regions]

      if (nrow(blast.best) == 0) {
        unlink(samp.dir, recursive = TRUE)
        print(paste0(samp, ": no contigs passed the read-depth filter."))
        writeLines(character(), out.fa)
        file.create(done.file)
        return(TRUE)
      }

      # Name contigs as region_contig_N (N = rank within region by bitscore desc)
      blast.best = blast.best[order(tName, -bitscore)]
      blast.best[, contig.idx := seq_len(.N), by = tName]
      blast.best[, new.name := paste0(tName, "_contig_", contig.idx)]

      # Subset and rename the DNAStringSet
      keep.names = blast.best$qName[blast.best$qName %in% names(ctgs)]
      all.contigs = ctgs[keep.names]
      names(all.contigs) = blast.best$new.name[match(keep.names, blast.best$qName)]

      # Write per-sample FASTA
      if (length(all.contigs) > 0) {
        Biostrings::writeXStringSet(all.contigs,
                                    filepath = paste0(output.directory, "/", samp, ".fa"))
        print(paste0(samp, ": assembled ", length(all.contigs), " contigs across ",
                     length(unique(blast.best$tName)), " regions."))
      } else {
        writeLines(character(), out.fa)
        print(paste0(samp, ": no contigs assembled."))
      }

      # Remove the per-sample working directory now that the FASTA is written
      unlink(samp.dir, recursive = TRUE)

      rm(all.contigs, blast.data, blast.best, ctgs)
      gc()
      file.create(done.file)
      TRUE

    }, error = function(e) {
      print(paste0("Error assembling ", sample.names[s], ": ", e$message))
      FALSE
    })
  }, mc.cores = workers)
  if (!all(vapply(sample.results, isTRUE, logical(1)))) {
    stop("Assembly failed for one or more samples. Correct the errors and rerun.")
  }

  # Clean up shared BLAST database
  unlink(Sys.glob(paste0(blast.db, ".*")))

  print("Shared region assembly complete.")

}#end function

#END SCRIPT
