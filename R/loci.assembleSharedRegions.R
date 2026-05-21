#' @title assembleSharedRegions
#'
#' @description Assembles per-sample contigs for each novel genomic region identified by
#' \code{discoverSharedRegions}. For each sample and each region in the shared-regions BED
#' file, reads are extracted from the sample's genome-mapped BAM, and SPAdes is run to
#' assemble contigs. Resulting contigs are named by region and combined into a single FASTA
#' per sample, written to \code{output.directory}. The output is structured to be compatible
#' with the downstream \code{filterHeterozygosity} → \code{annotateTargets} →
#' \code{alignTargets} pipeline (workflow 4), using the \code{novel_targets.fa} from
#' \code{discoverSharedRegions} as the target file.
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
#' @return Writes one FASTA file per sample to \code{output.directory}, where each contig
#' is named \code{region_contig_N} (e.g. \code{chr3_450000_450800_contig_1}). No value is
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
                                 samtools.path = NULL) {

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

  #Overwrite
  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
    } else {
      stop("Overwrite = FALSE and output directory exists. Either change to TRUE or overwrite manually.")
    }
  }
  dir.create(output.directory)

  #Load region BED file
  regions = data.table::fread(region.bed, sep = "\t", header = FALSE)
  data.table::setnames(regions, c("chrom", "start", "end", "sample_count"))
  region.names = paste0(regions$chrom, "_", regions$start, "_", regions$end)

  if (nrow(regions) == 0) { print("No regions found in novel_regions.bed."); return(NULL) }
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
      samp.dir = paste0(output.directory, "/", samp)
      dir.create(samp.dir, showWarnings = FALSE)

      all.contigs = Biostrings::DNAStringSet()

      for (r in 1:nrow(regions)) {

        reg.name = region.names[r]
        region.str = paste0(regions$chrom[r], ":", regions$start[r], "-", regions$end[r])
        tmp.bam = paste0(samp.dir, "/tmp_", reg.name, ".bam")
        tmp.fq  = paste0(samp.dir, "/tmp_", reg.name, ".fastq")
        spades.dir = paste0(samp.dir, "/spades_", reg.name)

        # Extract reads mapping to this region
        system(paste0(samtools.path, "samtools view -b ", bam,
                      " ", region.str, " -o ", tmp.bam),
               ignore.stdout = quiet, ignore.stderr = quiet)

        if (!file.exists(tmp.bam)) { next }

        n.reads = as.integer(trimws(
          system(paste0(samtools.path, "samtools view -c ", tmp.bam), intern = TRUE)))
        if (is.na(n.reads) || n.reads < min.reads.assemble) {
          system(paste0("rm -f ", tmp.bam))
          next
        }

        # Convert to FASTQ (single-end; mates outside region are not co-extracted)
        system(paste0(samtools.path, "samtools fastq ", tmp.bam, " > ", tmp.fq),
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

      # Write per-sample contig FASTA
      if (length(all.contigs) > 0) {
        Biostrings::writeXStringSet(all.contigs,
                                    filepath = paste0(output.directory, "/", samp, ".fa"))
        print(paste0(samp, ": assembled ", length(all.contigs), " contigs across ",
                     nrow(regions), " regions."))
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
