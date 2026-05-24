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
#' @param max.merge.distance maximum distance in bp between adjacent covered intervals within
#' a sample before they are merged into a single region. Default 500.
#'
#' @param min.mapping.quality minimum MAPQ score for a read to be counted. Filters
#' multi-mapping reads in repetitive regions. Default 20.
#'
#' @param threads number of CPU threads for parallel sample processing. Default 1.
#'
#' @param memory total RAM in GB available. Default 8.
#'
#' @param overwrite logical. If TRUE, deletes and recreates the output directory. Default FALSE.
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
#'   \item \code{sample-bams/} — per-sample genome-mapped BAM files
#'   \item \code{novel_regions.bed} — BED file of shared novel regions
#'   \item \code{novel_targets.fa} — genome sequences at shared regions (use as target file
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

  #Path normalization
  if (is.null(hisat2.path) == FALSE) {
    b.string = unlist(strsplit(hisat2.path, ""))
    if (b.string[length(b.string)] != "/") {
      hisat2.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { hisat2.path = "" }

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
  if (is.null(alignment.directory)) { print("alignment.directory not provided."); return(NULL) }
  if (is.null(read.directory)) { print("read.directory not provided."); return(NULL) }
  if (is.null(genome.file)) { print("genome.file not provided."); return(NULL) }
  if (is.null(output.directory)) { print("output.directory not provided."); return(NULL) }
  if (!file.exists(genome.file)) {
    print(paste0("genome.file not found: ", genome.file)); return(NULL)
  }

  # Verify external tools are reachable before doing any heavy work
  bt.check = Sys.which(paste0(bedtools.path, "bedtools"))
  if (bt.check == "") {
    stop("bedtools not found at '", bedtools.path, "bedtools'. ",
         "Run 'which bedtools' on the cluster and update bedtools.path in your config.")
  }

  # Create output directories; never wipe on resume — overwrite controls per-step redo
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(output.directory, "/sample-bams"), showWarnings = FALSE)
  dir.create(paste0(output.directory, "/covered-beds"), showWarnings = FALSE)

  #Gather alignments
  align.files = list.files(alignment.directory)
  if (length(align.files) == 0) { print("No alignment files found in alignment.directory."); return(NULL) }

  known.ref    = paste0(output.directory, "/known_loci_consensus.fa")
  known.index  = paste0(output.directory, "/known_loci_index")
  genome.index = paste0(output.directory, "/genome_index")

  ##################################################################################################
  ## Step 1: Build known-loci consensus reference
  ##################################################################################################
  if (overwrite == TRUE || !file.exists(known.ref) || file.size(known.ref) == 0) {
    print(paste0("Building known-loci consensus from ", length(align.files), " alignments..."))

    all.consensus = Biostrings::DNAStringSet()
    for (i in 1:length(align.files)) {
      locus.name = gsub("\\..*$", "", align.files[i])
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
    print("Known-loci consensus already exists — skipping Step 1.")
  }

  ##################################################################################################
  ## Step 2: Index references
  ##################################################################################################
  if (overwrite == TRUE || !file.exists(paste0(known.index, ".1.ht2")) || file.size(paste0(known.index, ".1.ht2")) == 0) {
    print("Indexing known-loci consensus and genome...")
    system(paste0(hisat2.path, "hisat2-build ", known.ref, " ", known.index),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(hisat2.path, "hisat2-build ", genome.file, " ", genome.index),
           ignore.stdout = quiet, ignore.stderr = quiet)
  } else {
    print("HISAT2 indices already exist — skipping Step 2.")
  }

  ##################################################################################################
  ## Step 3: Per-sample — map to known loci, extract unmapped, map to genome
  ##################################################################################################
  sample.names = list.dirs(read.directory, recursive = FALSE, full.names = FALSE)
  if (length(sample.names) == 0) {
    print("No sample directories found in read.directory."); return(NULL)
  }
  print(paste0("Processing ", length(sample.names), " samples..."))

  parallel::mclapply(sample.names, function(samp) {
    tryCatch({

      bed.out = paste0(output.directory, "/covered-beds/", samp, "_covered.bed")
      if (overwrite == FALSE && file.exists(bed.out) && file.size(bed.out) > 0) {
        print(paste0(samp, ": coverage BED already exists — skipping."))
        return(NULL)
      }

      read.dir = paste0(read.directory, "/", samp)
      if (!dir.exists(read.dir)) {
        print(paste0(samp, ": read directory not found. Skipping."))
        return(NULL)
      }

      read.files = list.files(read.dir, full.names = TRUE)
      r1 = read.files[grep("_R1_|_R1\\.f|-R1\\.f|READ1|_1\\.f|-1\\.f", read.files)]
      r2 = read.files[grep("_R2_|_R2\\.f|-R2\\.f|READ2|_2\\.f|-2\\.f", read.files)]
      if (length(r1) == 0 || length(r2) == 0) {
        print(paste0(samp, ": R1/R2 reads not found. Skipping."))
        return(NULL)
      }
      r1 = r1[1]; r2 = r2[1]

      tmp = paste0(output.directory, "/tmp_", samp)
      dir.create(tmp, showWarnings = FALSE)

      # Map to known loci, extract reads where BOTH mates are unmapped (-f 12)
      system(paste0(hisat2.path, "hisat2 -x ", known.index,
                    " -1 ", r1, " -2 ", r2,
                    " --threads 1 --no-spliced-alignment -S ", tmp, "/known.sam"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0(samtools.path, "samtools view -b -f 12 -F 256 ",
                    tmp, "/known.sam | ",
                    samtools.path, "samtools sort -n -o ", tmp, "/unmapped_sorted.bam"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0(samtools.path, "samtools fastq",
                    " -1 ", tmp, "/unmapped_R1.fastq",
                    " -2 ", tmp, "/unmapped_R2.fastq",
                    " -0 /dev/null -s /dev/null -n ",
                    tmp, "/unmapped_sorted.bam"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0("rm ", tmp, "/known.sam ", tmp, "/unmapped_sorted.bam"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      # Check that unmapped reads were produced
      if (!file.exists(paste0(tmp, "/unmapped_R1.fastq"))) {
        print(paste0(samp, ": no unmapped reads produced. Skipping."))
        system(paste0("rm -r ", tmp))
        return(NULL)
      }

      # Map unmapped reads to genome, keep mapped (-F 4), filter by MAPQ
      system(paste0(hisat2.path, "hisat2 -x ", genome.index,
                    " -1 ", tmp, "/unmapped_R1.fastq",
                    " -2 ", tmp, "/unmapped_R2.fastq",
                    " --threads 1 --no-spliced-alignment -S ", tmp, "/genome.sam"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0(samtools.path, "samtools view -b -F 4 -q ", min.mapping.quality,
                    " ", tmp, "/genome.sam | ",
                    samtools.path, "samtools sort -o ",
                    output.directory, "/sample-bams/", samp, ".bam"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0(samtools.path, "samtools index ",
                    output.directory, "/sample-bams/", samp, ".bam"),
             ignore.stdout = quiet, ignore.stderr = quiet)
      system(paste0("rm ", tmp, "/genome.sam"), ignore.stdout = quiet, ignore.stderr = quiet)

      # Get covered regions: depth >= min.coverage, merge nearby intervals
      bam.out = paste0(output.directory, "/sample-bams/", samp, ".bam")
      n.mapped = as.integer(trimws(
        system(paste0(samtools.path, "samtools view -c -F 4 ", bam.out), intern = TRUE)))

      if (is.na(n.mapped) || n.mapped == 0) {
        print(paste0(samp, ": no reads mapped to genome — covered BED will be empty."))
        file.create(bed.out)
        print(paste0("Finished mapping ", samp, " (0 genome-mapped reads)"))
      } else {
        # Pipe genomecov → awk filter → merge in one shot (no intermediate file on disk)
        # ignore.stderr = FALSE so bedtools errors are always visible
        system(paste0(bedtools.path, "bedtools genomecov -ibam ", bam.out, " -bg | ",
                      "awk '$4 >= ", min.coverage, "' | ",
                      bedtools.path, "bedtools merge -d ", max.merge.distance,
                      " > ", bed.out),
               ignore.stderr = FALSE)
        n.regions = as.integer(trimws(system(paste0("wc -l < ", bed.out), intern = TRUE)))
        print(paste0("Finished mapping ", samp, " (", n.mapped, " genome-mapped reads, ",
                     n.regions, " covered regions)"))
      }

      system(paste0("rm -r ", tmp))

    }, error = function(e) {
      print(paste0("Error processing ", samp, ": ", e$message))
      return(NULL)
    })
  }, mc.cores = threads)

  ##################################################################################################
  ## Step 4: Find genomic regions shared across >= min.samples samples
  ##################################################################################################
  shared.bed = paste0(output.directory, "/novel_regions.bed")
  novel.fa   = paste0(output.directory, "/novel_targets.fa")

  if (overwrite == TRUE || !file.exists(shared.bed) || file.size(shared.bed) == 0) {
    print("Finding shared genomic regions...")
    bed.files = list.files(paste0(output.directory, "/covered-beds"),
                           pattern = "_covered.bed$", full.names = TRUE)
    if (length(bed.files) == 0) {
      print("No per-sample coverage BEDs produced. Check mapping settings and read paths.")
      return(NULL)
    }

    # Concatenate per-sample BEDs with sample name in column 4 for count_distinct
    all.bed = paste0(output.directory, "/all_covered.bed")
    if (file.exists(all.bed)) { system(paste0("rm ", all.bed)) }
    non.empty.beds = bed.files[file.size(bed.files) > 0]
    print(paste0(length(non.empty.beds), " of ", length(bed.files),
                 " samples had covered genomic regions."))

    if (length(non.empty.beds) == 0) {
      print("All per-sample coverage BEDs are empty — no reads mapped to novel genomic regions.")
      print("Possible causes: reads fully explained by known loci; MAPQ threshold too strict;")
      print("or genome/read mismatch. Try lowering min.coverage or min.mapping.quality.")
      return(NULL)
    }

    for (i in 1:length(non.empty.beds)) {
      samp.label = gsub("_covered\\.bed$", "", basename(non.empty.beds[i]))
      system(paste0("awk -v s='", samp.label, "' 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,s}' ",
                    non.empty.beds[i], " >> ", all.bed))
    }

    # Sort, merge, count distinct samples per region, filter
    system(paste0("sort -k1,1 -k2,2n ", all.bed, " | ",
                  bedtools.path, "bedtools merge -i stdin -c 4 -o count_distinct | ",
                  "awk -v ms=", min.samples, " -v ml=", min.region.length,
                  " '$4 >= ms && ($3-$2) >= ml' > ", shared.bed))
    system(paste0("rm ", all.bed))

    n.regions = as.integer(trimws(system(paste0("wc -l < ", shared.bed), intern = TRUE)))
    if (is.na(n.regions) || n.regions == 0) {
      print(paste0("No shared novel regions found in >= ", min.samples, " samples."))
      print("Try lowering min.samples, min.coverage, or min.region.length.")
      return(NULL)
    }
    print(paste0("Found ", n.regions, " novel regions covered in >= ", min.samples, " samples."))
  } else {
    print("Shared regions BED already exists — skipping Step 4.")
  }

  ##################################################################################################
  ## Step 5: Extract genome sequences at shared regions → novel_targets.fa
  ##################################################################################################
  if (overwrite == TRUE || !file.exists(novel.fa) || file.size(novel.fa) == 0) {
    raw.fa = paste0(output.directory, "/novel_targets_raw.fa")
    system(paste0(bedtools.path, "bedtools getfasta -fi ", genome.file,
                  " -bed ", shared.bed, " -fo ", raw.fa),
           ignore.stdout = quiet, ignore.stderr = quiet)

    raw.seqs = Biostrings::readDNAStringSet(raw.fa)
    # Rename chr:start-end → chr_start_end (phylip-safe, no colons)
    # Sanitize sequence names: bedtools getfasta produces "chrom:start-end".
    # Replace ":" and "-" first (coordinate delimiters), then replace any remaining
    # shell-unsafe character (;  =  space  etc.) with "_" so names are safe to use
    # in file paths and shell commands.
    names(raw.seqs) = gsub("[^A-Za-z0-9_.]", "_", gsub("[:-]", "_", names(raw.seqs)))
    Biostrings::writeXStringSet(raw.seqs, filepath = novel.fa)
    system(paste0("rm ", raw.fa))
    rm(raw.seqs)

    # Cleanup indices and temp consensus (only after first successful run)
    system(paste0("rm -f ", known.ref, " ", known.index, ".* ", genome.index, ".*"))
  } else {
    print("Novel targets FASTA already exists — skipping Step 5.")
  }

  print(paste0("Novel target sequences written to: ", novel.fa))
  print(paste0("Per-sample BAM files written to:   ", output.directory, "/sample-bams/"))
  print(paste0("Shared regions BED:                ", shared.bed))

}#end function

#END SCRIPT
