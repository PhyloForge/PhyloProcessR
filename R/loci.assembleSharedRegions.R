#' @title assembleSharedRegions
#'
#' @description Assembles per-sample contigs for each novel genomic region identified by
#' \code{discoverSharedRegions}. For each sample, reads overlapping any novel region are
#' extracted from the genome-mapped BAM in a single pass, then a read-count filter is applied
#' across all regions simultaneously (via \code{bedtools coverage}) so that SPAdes is only
#' invoked for regions that actually have sufficient coverage. Resulting contigs are named by
#' region and combined into a single FASTA per sample, written to \code{output.directory}. The
#' output is structured to be compatible with the downstream
#' \code{filterHeterozygosity} → \code{collectNovelContigs} → \code{alignTargets} pipeline
#' (workflow X4).
#'
#' @param discover.directory path to the output directory from \code{discoverSharedRegions}.
#' Must contain \code{sample-bams/} and \code{novel_regions.bed}.
#'
#' @param output.directory path to write per-sample contig FASTA files
#' (e.g. \code{data-analysis/contigs/9_genome-contigs}).
#'
#' @param min.reads.assemble minimum number of reads mapping to a region for SPAdes to be
#' run. Regions with fewer reads are skipped. Default 5.
#'
#' @param kmer.values integer vector of k-mer sizes passed to SPAdes. Default
#' \code{c(33, 55, 77, 99, 127)}.
#'
#' @param threads number of CPU threads for parallel sample processing. Default 1.
#'
#' @param memory total RAM in GB available. Default 8.
#'
#' @param overwrite logical. If TRUE, deletes and recreates the output directory. Default FALSE.
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
#' @return Writes one FASTA file per sample to \code{output.directory}, where each contig
#' is named \code{region_contig_N} (e.g. \code{chr3_450000_450800_contig_1}). Intermediate
#' per-sample working directories (temp BAMs, FASTQs, SPAdes output) are written under
#' \code{discover.directory} and cleaned up after each sample completes. No value is
#' returned to R.
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
                                 bedtools.path = NULL) {

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

  #Path normalization
  if (is.null(spades.path) == FALSE) {
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { spades.path = "" }

  if (is.null(samtools.path) == FALSE) {
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { samtools.path = "" }

  if (is.null(bedtools.path) == FALSE) {
    b.string = unlist(strsplit(bedtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      bedtools.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { bedtools.path = "" }

  #Input checks
  if (is.null(discover.directory)) { print("discover.directory not provided."); return(NULL) }
  if (is.null(output.directory)) { print("output.directory not provided."); return(NULL) }

  region.bed = paste0(discover.directory, "/novel_regions.bed")
  bam.dir = paste0(discover.directory, "/sample-bams")
  if (!file.exists(region.bed)) {
    print(paste0("novel_regions.bed not found in: ", discover.directory)); return(NULL)
  }
  if (!dir.exists(bam.dir)) {
    print(paste0("sample-bams/ directory not found in: ", discover.directory)); return(NULL)
  }

  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)

  #Load region BED file
  regions = data.table::fread(region.bed, sep = "\t", header = FALSE)
  if (is.null(regions) || nrow(regions) == 0) {
    print("novel_regions.bed is empty or missing — no shared regions were found.")
    print("Check that bedtools ran successfully in discoverSharedRegions (correct bedtools.path?)")
    return(NULL)
  }
  data.table::setnames(regions, c("chrom", "start", "end", "sample_count"))
  region.names = paste0(regions$chrom, "_", regions$start, "_", regions$end)
  print(paste0("Assembling contigs for ", nrow(regions), " regions..."))

  #Get per-sample BAM files
  bam.files = list.files(bam.dir, pattern = "\\.bam$", full.names = TRUE)
  bam.files = bam.files[!grepl("\\.bai$", bam.files)]
  sample.names = gsub("\\.bam$", "", basename(bam.files))

  if (length(bam.files) == 0) { print("No BAM files found in sample-bams/."); return(NULL) }
  print(paste0("Assembling ", length(sample.names), " samples..."))

  kmer.str = paste0(kmer.values, collapse = ",")
  mem.per = max(floor(memory / threads), 4)

  ##################################################################################################
  ## Per-sample assembly (parallel)
  ##################################################################################################
  parallel::mclapply(seq_along(sample.names), function(s) {
    tryCatch({

      samp = sample.names[s]
      bam = bam.files[s]

      if (overwrite == FALSE && file.exists(paste0(output.directory, "/", samp, ".fa"))) {
        print(paste0(samp, ": contig FASTA already exists — skipping."))
        return(NULL)
      }

      # Intermediate work dir goes inside discover.directory to keep output.directory clean
      samp.dir = paste0(discover.directory, "/", samp)
      dir.create(samp.dir, showWarnings = FALSE)

      ##########################################################################
      # Step A: Extract all reads overlapping any novel region in one BAM pass.
      # Working from this small filtered BAM instead of the full genome BAM
      # makes every subsequent operation orders of magnitude faster.
      ##########################################################################
      novel.bam = paste0(samp.dir, "/novel_reads.bam")
      ret0 = system(paste0(samtools.path, "samtools view -b -L ", region.bed,
                           " -o ", novel.bam, " ", bam),
                    ignore.stdout = quiet, ignore.stderr = quiet)
      if (ret0 != 0 || !file.exists(novel.bam) || file.size(novel.bam) == 0) {
        print(paste0(samp, ": failed to extract novel-region reads. Skipping."))
        return(NULL)
      }
      system(paste0(samtools.path, "samtools index ", novel.bam),
             ignore.stdout = quiet, ignore.stderr = quiet)

      n.novel.reads = as.integer(trimws(
        system(paste0(samtools.path, "samtools view -c ", novel.bam), intern = TRUE)))
      print(paste0(samp, ": ", n.novel.reads, " reads mapped to novel regions."))

      if (is.na(n.novel.reads) || n.novel.reads == 0) {
        system(paste0("rm -f ", novel.bam, " ", novel.bam, ".bai"))
        print(paste0(samp, ": no reads in novel regions — skipping."))
        return(NULL)
      }

      ##########################################################################
      # Step B: Count reads per region in one bedtools pass.
      # bedtools coverage -counts appends one column: number of reads overlapping
      # each BED interval. This replaces 26k+ individual samtools view calls.
      ##########################################################################
      cov.file = paste0(samp.dir, "/region_readcounts.tsv")
      system(paste0(bedtools.path, "bedtools coverage -a ", region.bed,
                    " -b ", novel.bam, " -counts > ", cov.file),
             ignore.stdout = quiet, ignore.stderr = quiet)

      region.counts = data.table::fread(cov.file, header = FALSE)
      reads.per.region = region.counts[[ncol(region.counts)]]
      system(paste0("rm -f ", cov.file))

      active.r = which(reads.per.region >= min.reads.assemble)
      print(paste0(samp, ": ", length(active.r), " of ", nrow(regions),
                   " regions have >= ", min.reads.assemble, " reads — running SPAdes."))

      if (length(active.r) == 0) {
        system(paste0("rm -f ", novel.bam, " ", novel.bam, ".bai"))
        print(paste0(samp, ": no regions with sufficient reads."))
        return(NULL)
      }

      ##########################################################################
      # Step C: Per-region SPAdes assembly (only active regions)
      ##########################################################################
      all.contigs = Biostrings::DNAStringSet()

      for (r in active.r) {

        reg.name = region.names[r]
        # BED is 0-based half-open; samtools view uses 1-based closed — add 1 to start
        region.str = paste0(regions$chrom[r], ":", regions$start[r] + 1L, "-", regions$end[r])
        tmp.bam    = paste0(samp.dir, "/tmp_", reg.name, ".bam")
        tmp.fq     = paste0(samp.dir, "/tmp_", reg.name, ".fastq")
        spades.dir = paste0(samp.dir, "/spades_", reg.name)

        # Extract reads for this region from the pre-filtered novel BAM
        ret = system(paste0(samtools.path, "samtools view -b -o ", tmp.bam,
                            " ", novel.bam, " ", region.str),
                     ignore.stdout = quiet, ignore.stderr = quiet)
        if (ret != 0 || !file.exists(tmp.bam)) { next }

        # Convert to FASTQ
        system(paste0(samtools.path, "samtools fastq -o ", tmp.fq, " ", tmp.bam),
               ignore.stdout = quiet, ignore.stderr = quiet)

        # SPAdes assembly
        system(paste0(spades.path, "spades.py -s ", tmp.fq,
                      " -k ", kmer.str,
                      " -o ", spades.dir,
                      " -m ", mem.per,
                      " --threads 1 --careful"),
               ignore.stdout = quiet, ignore.stderr = quiet)

        # Load contigs if assembly produced output
        contig.fa = paste0(spades.dir, "/contigs.fasta")
        if (file.exists(contig.fa)) {
          ctgs = Biostrings::readDNAStringSet(contig.fa)
          if (length(ctgs) > 0) {
            names(ctgs) = paste0(reg.name, "_contig_", seq_along(ctgs))
            all.contigs = append(all.contigs, ctgs)
          }
        }

        system(paste0("rm -rf ", spades.dir, " ", tmp.bam, " ", tmp.fq))

      }#end region loop

      # Clean up the pre-filtered novel BAM
      system(paste0("rm -f ", novel.bam, " ", novel.bam, ".bai"))

      # Write per-sample contig FASTA
      if (length(all.contigs) > 0) {
        Biostrings::writeXStringSet(all.contigs,
                                    filepath = paste0(output.directory, "/", samp, ".fa"))
        print(paste0(samp, ": assembled ", length(all.contigs), " contigs across ",
                     length(active.r), " regions."))
      } else {
        print(paste0(samp, ": no contigs assembled."))
      }

      rm(all.contigs)
      gc()

    }, error = function(e) {
      print(paste0("Error assembling ", sample.names[s], ": ", e$message))
      return(NULL)
    })
  }, mc.cores = threads)

  print("Shared region assembly complete.")

}#end function

#END SCRIPT
