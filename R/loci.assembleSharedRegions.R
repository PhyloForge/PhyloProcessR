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
#' \code{filterHeterozygosity} → \code{collectNovelContigs} → \code{alignTargets} pipeline
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

  # Path normalisation
  norm.path = function(p) {
    if (is.null(p)) return("")
    b = unlist(strsplit(p, ""))
    if (b[length(b)] != "/") p = paste0(p, "/")
    p
  }
  spades.path   = norm.path(spades.path)
  samtools.path = norm.path(samtools.path)
  bedtools.path = norm.path(bedtools.path)
  blast.path    = norm.path(blast.path)

  # Input checks
  if (is.null(discover.directory)) { print("discover.directory not provided."); return(NULL) }
  if (is.null(output.directory))   { print("output.directory not provided.");   return(NULL) }

  region.bed  = paste0(discover.directory, "/novel_regions.bed")
  novel.fa    = paste0(discover.directory, "/novel_targets.fa")
  bam.dir     = paste0(discover.directory, "/sample-bams")

  if (!file.exists(region.bed)) {
    print(paste0("novel_regions.bed not found in: ", discover.directory)); return(NULL)
  }
  if (!file.exists(novel.fa)) {
    print(paste0("novel_targets.fa not found in: ", discover.directory)); return(NULL)
  }
  if (!dir.exists(bam.dir)) {
    print(paste0("sample-bams/ not found in: ", discover.directory)); return(NULL)
  }

  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)

  # Load region BED
  regions = data.table::fread(region.bed, sep = "\t", header = FALSE)
  if (is.null(regions) || nrow(regions) == 0) {
    print("novel_regions.bed is empty — no shared regions were found.")
    return(NULL)
  }
  data.table::setnames(regions, c("chrom", "start", "end", "sample_count"))
  # Sanitize region names the same way discoverSharedRegions sanitizes novel_targets.fa
  # sequence names — must match exactly so the BLAST tName lookup works.
  # Scaffold names from some genome assemblies contain shell-unsafe characters
  # (e.g. "ScWFrIx_100724;HRSCAF=196187") that break file paths and MAFFT commands.
  region.names = gsub("[^A-Za-z0-9_.]", "_",
                      paste0(regions$chrom, "_", regions$start, "_", regions$end))
  print(paste0("Assembling contigs for ", nrow(regions), " regions..."))

  # Get per-sample BAM files
  bam.files = list.files(bam.dir, pattern = "\\.bam$", full.names = TRUE)
  bam.files = bam.files[!grepl("\\.bai$", bam.files)]
  sample.names = gsub("\\.bam$", "", basename(bam.files))

  if (length(bam.files) == 0) { print("No BAM files found in sample-bams/."); return(NULL) }
  print(paste0("Assembling ", length(sample.names), " samples..."))

  ##################################################################################################
  ## Build a shared BLAST database from novel_targets.fa (once, before parallel loop)
  ##################################################################################################
  blast.db = paste0(discover.directory, "/novel_targets_blast_db")
  system(paste0(blast.path, "makeblastdb -in ", novel.fa,
                " -dbtype nucl -out ", blast.db),
         ignore.stdout = quiet, ignore.stderr = quiet)

  kmer.str  = paste0(kmer.values, collapse = ",")
  mem.per   = max(floor(memory / threads), 4)
  blast.headers = c("qName", "tName", "pident", "length", "mismatch", "gapopen",
                    "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore")

  ##################################################################################################
  ## Per-sample assembly (parallel)
  ##################################################################################################
  parallel::mclapply(seq_along(sample.names), function(s) {
    tryCatch({

      samp = sample.names[s]
      bam  = bam.files[s]

      if (overwrite == FALSE && file.exists(paste0(output.directory, "/", samp, ".fa"))) {
        print(paste0(samp, ": contig FASTA already exists — skipping."))
        return(NULL)
      }

      samp.dir = paste0(discover.directory, "/", samp)
      dir.create(samp.dir, showWarnings = FALSE)

      ##########################################################################
      # Step A: Extract all reads overlapping any novel region in one BAM pass
      ##########################################################################
      novel.bam = paste0(samp.dir, "/novel_reads.bam")
      ret0 = system(paste0(samtools.path, "samtools view -b -L ", region.bed,
                           " -o ", novel.bam, " ", bam),
                    ignore.stdout = quiet, ignore.stderr = quiet)
      if (ret0 != 0 || !file.exists(novel.bam) || file.size(novel.bam) == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": failed to extract novel-region reads. Skipping."))
        return(NULL)
      }
      system(paste0(samtools.path, "samtools index ", novel.bam),
             ignore.stdout = quiet, ignore.stderr = quiet)

      n.novel = as.integer(trimws(
        system(paste0(samtools.path, "samtools view -c ", novel.bam), intern = TRUE)))
      print(paste0(samp, ": ", n.novel, " reads mapped to novel regions."))

      if (is.na(n.novel) || n.novel == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": no reads in novel regions — skipping."))
        return(NULL)
      }

      ##########################################################################
      # Step B: Count reads per region in one bedtools pass to determine which
      # regions have sufficient coverage (used to filter contigs after BLAST)
      ##########################################################################
      cov.out = system(paste0(bedtools.path, "bedtools coverage -a ", region.bed,
                              " -b ", novel.bam, " -counts"),
                       intern = TRUE, ignore.stderr = quiet)

      if (length(cov.out) == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": bedtools coverage produced no output — skipping."))
        return(NULL)
      }

      region.counts    = data.table::fread(text = paste(cov.out, collapse = "\n"), header = FALSE)
      if (ncol(region.counts) == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": bedtools coverage output could not be parsed — skipping."))
        return(NULL)
      }
      reads.per.region = region.counts[[ncol(region.counts)]]
      active.regions   = region.names[reads.per.region >= min.reads.assemble]

      print(paste0(samp, ": ", length(active.regions), " of ", nrow(regions),
                   " regions have >= ", min.reads.assemble, " reads."))

      if (length(active.regions) == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": no regions with sufficient reads — skipping."))
        return(NULL)
      }

      ##########################################################################
      # Step C: Convert all novel reads to FASTQ and run ONE SPAdes assembly
      ##########################################################################
      novel.fq   = paste0(samp.dir, "/novel_reads.fastq")
      spades.dir = paste0(samp.dir, "/spades_all")

      system(paste0(samtools.path, "samtools fastq -o ", novel.fq, " ", novel.bam),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0("rm -f ", novel.bam, " ", novel.bam, ".bai"))

      system(paste0(spades.path, "spades.py -s ", novel.fq,
                    " -k ", kmer.str,
                    " -o ", spades.dir,
                    " -m ", mem.per,
                    " --threads 1 --careful"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0("rm -f ", novel.fq))

      contig.fa = paste0(spades.dir, "/contigs.fasta")
      if (!file.exists(contig.fa) || file.size(contig.fa) == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": SPAdes produced no contigs."))
        return(NULL)
      }

      ##########################################################################
      # Step D: BLAST contigs against novel_targets.fa to assign region names.
      # tName values in the output are chr_start_end, matching region.names.
      ##########################################################################
      blast.out = paste0(samp.dir, "/contig_blast.txt")
      system(paste0(blast.path, "blastn -task blastn -db ", blast.db,
                    " -query ", contig.fa,
                    " -out ", blast.out,
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen",
                    " qstart qend sstart send evalue bitscore\"",
                    " -num_threads 1 -evalue 0.001"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      ctgs = Biostrings::readDNAStringSet(contig.fa)
      system(paste0("rm -rf ", spades.dir))

      if (!file.exists(blast.out) || file.size(blast.out) == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": no BLAST hits for assembled contigs."))
        return(NULL)
      }

      blast.data = data.table::fread(blast.out, header = FALSE)
      system(paste0("rm -f ", blast.out))

      if (nrow(blast.data) == 0) {
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": no BLAST hits for assembled contigs."))
        return(NULL)
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
        system(paste0("rm -rf ", samp.dir))
        print(paste0(samp, ": no contigs passed the read-depth filter."))
        return(NULL)
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
        print(paste0(samp, ": no contigs assembled."))
      }

      # Remove the per-sample working directory now that the FASTA is written
      system(paste0("rm -rf ", samp.dir))

      rm(all.contigs, blast.data, blast.best, ctgs)
      gc()

    }, error = function(e) {
      print(paste0("Error assembling ", sample.names[s], ": ", e$message))
      return(NULL)
    })
  }, mc.cores = threads)

  # Clean up shared BLAST database
  system(paste0("rm -f ", blast.db, "*"))

  print("Shared region assembly complete.")

}#end function

#END SCRIPT
