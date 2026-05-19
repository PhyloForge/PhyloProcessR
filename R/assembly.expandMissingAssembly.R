#' @title expandMissingAssembly
#'
#' @description Attempts to recover target loci that are missing or poorly
#'   assembled by (1) BLASTing each sample's existing assembly against the
#'   reference to find matched loci, (2) computing a cross-sample consensus of
#'   the best contig per target locus, (3) mapping each sample's raw reads to
#'   its missing targets with HISAT2, (4) extracting mapped reads with samtools,
#'   (5) merging paired-end reads with fastp, and (6) assembling de novo with
#'   SPAdes. Newly assembled contigs are BLASTed back against the missing
#'   targets; those that pass quality filters are combined with the original
#'   matching contigs into a final expanded assembly. Phase 1 (BLAST matching)
#'   runs in parallel across samples; Phase 2 (mapping + assembly) runs
#'   sequentially per sample but passes full thread counts to each tool.
#'
#' @param assembly.directory path to a directory of per-sample assembly FASTA
#'   files (one file per sample, named \code{<sample>.fa}).
#'
#' @param read.directory path to the top-level processed reads directory (e.g.
#'   \code{"processed-reads"}). The actual reads used for mapping are taken from
#'   the subdirectory specified by \code{mapping.reads}.
#'
#' @param mapping.reads name of the subdirectory within \code{read.directory}
#'   that contains the per-sample read folders to use for Phase 2 read mapping.
#'   Should be paired (non-merged) reads. Common choices: \code{"decontaminated-reads"}
#'   (default), \code{"cleaned-reads"}, \code{"pe-merged-reads"}. Using
#'   pre-merged reads is not recommended for mapping as HISAT2 expects paired
#'   input.
#'
#' @param reference path to the FASTA file of reference/target sequences used
#'   for BLAST matching and HISAT2 mapping.
#'
#' @param output.directory path to the directory where per-sample working files
#'   and the final expanded assemblies will be written. Final assemblies are
#'   saved to \code{output.directory/expanded-assemblies/}. Default:
#'   \code{"expand-missing-assembly"}.
#'
#' @param min.match.percent minimum BLAST percent identity required to accept a
#'   contig-to-target match. Default: \code{60}.
#'
#' @param min.match.length minimum BLAST alignment length (bp) required to
#'   accept a match. Default: \code{100}.
#'
#' @param min.match.coverage minimum percentage of the target length that a
#'   BLAST match must cover to be accepted. Default: \code{35}.
#'
#' @param phase2.reference what to use as the HISAT2 mapping scaffold for loci
#'   that are missing from a given sample but present in at least one other
#'   sample. \code{"contig"} (default) uses the best assembled contig for that
#'   target from across all other samples — a closer match to the true sequence,
#'   giving better read recovery. \code{"reference"} uses the original probe/
#'   bait sequence from the reference file — useful when cross-sample contigs
#'   are unavailable or of low quality.
#'
#' @param recover.all.missing logical; if \code{TRUE}, a second recovery pass
#'   is performed for loci that are absent from every sample's assembly (i.e.
#'   not present in \code{final.save} at all). These loci are always mapped
#'   against the original reference sequences since no cross-sample contig
#'   exists. This pass can be time-consuming for large datasets. Default:
#'   \code{FALSE}.
#'
#' @param memory RAM in GB to pass to SPAdes and fastp. Default: \code{8}.
#'
#' @param threads number of CPU threads passed to BLAST (Phase 1 parallelism),
#'   and to HISAT2, samtools, and SPAdes within each Phase 2 sample. Default:
#'   \code{1}.
#'
#' @param spades.path path to the directory containing \code{spades.py}. If
#'   \code{NULL}, expected on the system PATH. Default: \code{NULL}.
#'
#' @param hisat2.path path to the directory containing \code{hisat2} and
#'   \code{hisat2-build}. If \code{NULL}, expected on the system PATH. Default:
#'   \code{NULL}.
#'
#' @param samtools.path path to the directory containing \code{samtools}. If
#'   \code{NULL}, expected on the system PATH. Default: \code{NULL}.
#'
#' @param fastp.path path to the directory containing \code{fastp}. If
#'   \code{NULL}, expected on the system PATH. Default: \code{NULL}.
#'
#' @param blast.path path to the directory containing the BLAST executables
#'   (\code{makeblastdb}, \code{blastn}). If \code{NULL}, expected on the
#'   system PATH. Default: \code{NULL}.
#'
#' @param overwrite logical; if \code{TRUE} already-completed samples in
#'   \code{output.directory/expanded-assemblies/} are overwritten. Default:
#'   \code{FALSE}.
#'
#' @param quiet logical; if \code{TRUE} tool stdout/stderr is suppressed.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Expanded assemblies are saved as
#'   \code{output.directory/expanded-assemblies/<sample>.fa}.
#'
#' @export

expandMissingAssembly = function(assembly.directory = NULL,
                                 read.directory = NULL,
                                 mapping.reads = "decontaminated-reads",
                                 reference = NULL,
                                 output.directory = "expand-missing-assembly",
                                 min.match.percent = 60,
                                 min.match.length = 100,
                                 min.match.coverage = 35,
                                 phase2.reference = c("contig", "reference"),
                                 recover.all.missing = FALSE,
                                 memory = 8,
                                 threads = 1,
                                 spades.path = NULL,
                                 hisat2.path = NULL,
                                 samtools.path = NULL,
                                 fastp.path = NULL,
                                 blast.path = NULL,
                                 overwrite = FALSE,
                                 quiet = TRUE) {

  # # Debug
  # library(PhyloProcessR)
  # setwd("/Volumes/LaCie/Mantellidae")
  # assembly.directory = "/Volumes/LaCie/Mantellidae/draft-assemblies"
  # output.directory   = "expand-missing-assembly"
  # reference = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
  # read.directory  = "/Volumes/LaCie/Mantellidae/processed-reads"
  # mapping.reads   = "decontaminated-reads"
  # spades.path    = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # samtools.path  = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # hisat2.path    = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # blast.path     = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # fastp.path     = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # quiet = TRUE
  # overwrite = FALSE
  # threads = 8
  # memory = 24
  # min.match.percent  = 60
  # min.match.length   = 100
  # min.match.coverage = 35
  # phase2.reference   = "contig"
  # recover.all.missing = FALSE

  phase2.reference = match.arg(phase2.reference)

  # Normalize tool paths
  norm.path = function(p) {
    if (is.null(p)) return("")
    b = unlist(strsplit(p, ""))
    if (b[length(b)] != "/") p = paste0(p, "/")
    p
  }
  blast.path    = norm.path(blast.path)
  hisat2.path   = norm.path(hisat2.path)
  samtools.path = norm.path(samtools.path)
  fastp.path    = norm.path(fastp.path)
  spades.path   = norm.path(spades.path)

  # Parameter checks
  if (is.null(assembly.directory)) stop("assembly.directory is required.")
  if (is.null(read.directory))     stop("read.directory is required.")
  if (is.null(reference))          stop("reference is required.")
  if (!dir.exists(assembly.directory)) stop("assembly.directory not found.")
  if (!dir.exists(read.directory))     stop("read.directory not found.")
  if (!file.exists(reference))         stop("reference file not found.")

  actual.read.dir = paste0(read.directory, "/", mapping.reads)
  if (!dir.exists(actual.read.dir)) {
    stop("mapping.reads subdirectory not found: ", actual.read.dir,
         "\n  Check that mapping.reads matches a subdirectory of read.directory.")
  }

  # Output directories
  expanded.dir = paste0(output.directory, "/expanded-assemblies")
  if (!dir.exists(output.directory)) dir.create(output.directory, recursive = TRUE)
  if (!dir.exists(expanded.dir))     dir.create(expanded.dir)

  file.names = list.files(assembly.directory)
  read.sets  = list.files(actual.read.dir)

  if (length(file.names) == 0) stop("No assembly files found in assembly.directory.")

  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  ###############################################################################
  ## Build BLAST reference DB once (major speedup — was rebuilt per sample)
  ###############################################################################
  db.dir = paste0(output.directory, "/blast_ref_db")
  if (!dir.exists(db.dir)) dir.create(db.dir)
  system(paste0(blast.path, "makeblastdb -in ", reference,
                " -parse_seqids -dbtype nucl -out ", db.dir, "/nucl-blast_db"),
         ignore.stdout = quiet, ignore.stderr = quiet)

  ###############################################################################
  ## Phase 1: BLAST each assembly against reference — parallelized
  ###############################################################################
  parallel::mclapply(seq_along(file.names), function(i) {
  tryCatch({

    sample      = gsub("\\..*$", "", file.names[i])
    species.dir = paste0(output.directory, "/", sample)
    if (!dir.exists(species.dir)) dir.create(species.dir)

    # Skip if already done
    if (overwrite == FALSE &&
        file.exists(paste0(species.dir, "/filtered-blast-match.txt"))) {
      print(paste0(sample, " Phase 1 already done, skipping."))
      return(NULL)
    }

    # BLAST assembly against shared reference DB
    blast.out = paste0(species.dir, "/target-blast-match.txt")
    system(paste0(
      blast.path, "blastn -task dc-megablast -db ", db.dir, "/nucl-blast_db -evalue 0.001",
      " -query ", assembly.directory, "/", file.names[i],
      " -out ", blast.out,
      " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend",
      " sstart send evalue bitscore qlen slen gaps\"",
      " -num_threads 1"
    ), ignore.stdout = quiet, ignore.stderr = quiet)

    if (!file.exists(blast.out) || file.size(blast.out) == 0) {
      print(paste0(sample, ": no BLAST matches to reference. Skipping."))
      return(NULL)
    }

    match.data = data.table::fread(blast.out, sep = "\t", header = FALSE,
                                   stringsAsFactors = FALSE)
    data.table::setnames(match.data, headers)
    system(paste0("rm -f ", blast.out))

    # Filter by identity, length, and coverage
    filt.data = match.data[match.data$matches > min.match.length, ]
    filt.data = filt.data[filt.data$pident >= min.match.percent, ]
    filt.data = filt.data[filt.data$matches >= ((min.match.coverage / 100) * filt.data$tLen), ]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, ": no matches passed filters. Skipping."))
      return(NULL)
    }

    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    sample.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", file.names[i]))
    match.contigs  = sample.contigs[names(sample.contigs) %in% filt.data$qName]

    Biostrings::writeXStringSet(match.contigs,
                                paste0(species.dir, "/", sample, "_matching-contigs.fa"))
    write.table(filt.data,
                file = paste0(species.dir, "/filtered-blast-match.txt"),
                row.names = FALSE, quote = FALSE, sep = "\t")

    print(paste0(sample, " Phase 1 complete: ",
                 length(unique(filt.data$tName)), " targets matched."))

  }, error = function(e) {
    warning(file.names[i], " Phase 1 failed: ", conditionMessage(e))
  })
  }, mc.cores = threads) # end Phase 1 mclapply

  ###############################################################################
  ## Cross-sample deduplication: best contig per target across all samples
  ###############################################################################
  all.matches = Biostrings::DNAStringSet()

  for (i in seq_along(file.names)) {
    sample      = gsub("\\..*$", "", file.names[i])
    species.dir = paste0(output.directory, "/", sample)
    match.file  = paste0(species.dir, "/", sample, "_matching-contigs.fa")
    blast.file  = paste0(species.dir, "/filtered-blast-match.txt")

    if (!file.exists(match.file) || !file.exists(blast.file)) next

    match.contigs = Biostrings::readDNAStringSet(match.file)
    filt.data     = read.table(blast.file, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)

    for (j in seq_len(nrow(filt.data))) {
      temp.match = match.contigs[names(match.contigs) %in% filt.data$qName[j]]
      if (length(temp.match) == 0) next
      names(temp.match) = filt.data$tName[j]
      all.matches = append(all.matches, temp.match)
    }
  }

  # Keep longest contig per target name
  dup.names  = unique(names(all.matches)[duplicated(names(all.matches))])
  final.save = all.matches[!names(all.matches) %in% dup.names]

  for (j in seq_along(dup.names)) {
    dup.contigs = all.matches[names(all.matches) %in% dup.names[j]]
    best        = dup.contigs[Biostrings::width(dup.contigs) == max(Biostrings::width(dup.contigs))][1]
    final.save  = append(final.save, best)
  }

  Biostrings::writeXStringSet(final.save,
                              paste0(output.directory, "/unique_matches.fa"))

  if (length(final.save) == 0 && recover.all.missing == FALSE) {
    message("No cross-sample matches found. Cannot proceed to Phase 2.")
    system(paste0("rm -rf ", db.dir))
    return(invisible(NULL))
  }

  # Load reference sequences if needed for phase2.reference="reference" or
  # recover.all.missing=TRUE (lazy — only pay the I/O cost when required)
  if (phase2.reference == "reference" || recover.all.missing == TRUE) {
    reference.seqs = Biostrings::readDNAStringSet(reference)
  }

  ###############################################################################
  ## Phase 2: Map reads to missing targets, assemble, expand — sequential so
  ##           each tool (HISAT2, SPAdes, samtools) gets the full thread budget
  ###############################################################################
  for (i in seq_along(file.names)) {
  tryCatch({

    sample      = gsub("\\..*$", "", file.names[i])
    species.dir = paste0(output.directory, "/", sample)
    out.file    = paste0(expanded.dir, "/", sample, ".fa")

    # Skip if already done
    if (overwrite == FALSE && file.exists(out.file)) {
      print(paste0(sample, " already finished, skipping."))
      next
    }

    blast.file = paste0(species.dir, "/filtered-blast-match.txt")
    if (!file.exists(blast.file)) {
      print(paste0(sample, ": no Phase 1 results found, skipping."))
      next
    }

    found.data   = read.table(blast.file, sep = "\t", header = TRUE,
                              stringsAsFactors = FALSE)
    found.targets = unique(found.data$tName)

    # Build the mapping reference for Phase 2 based on user settings
    if (phase2.reference == "contig") {
      # Default: use the best cross-sample contig for each missing target
      mapping.ref = final.save[!names(final.save) %in% found.targets]

      # Optionally also recover loci absent from ALL samples using original ref
      if (recover.all.missing == TRUE) {
        universal.names = names(reference.seqs)[!names(reference.seqs) %in% names(final.save)]
        universal.missing = reference.seqs[names(reference.seqs) %in% universal.names]
        if (length(universal.missing) > 0) {
          mapping.ref = append(mapping.ref, universal.missing)
          print(paste0(sample, ": adding ", length(universal.missing),
                       " universally missing targets from original reference."))
        }
      }
    } else {
      # phase2.reference == "reference": use original reference sequences for
      # all missing targets (cross-sample and universal covered together)
      missing.names = names(reference.seqs)[!names(reference.seqs) %in% found.targets]
      mapping.ref   = reference.seqs[names(reference.seqs) %in% missing.names]
    }

    # If all targets were already found, just save originals
    if (length(mapping.ref) == 0) {
      print(paste0(sample, ": all targets already assembled — saving originals."))
      found.contigs = Biostrings::readDNAStringSet(
        paste0(species.dir, "/", sample, "_matching-contigs.fa"))
      Biostrings::writeXStringSet(found.contigs, out.file)
      next
    }

    missing.ref = paste0(species.dir, "/missing_ref.fa")
    Biostrings::writeXStringSet(mapping.ref, missing.ref)

    # Gather reads for this sample
    input.reads = read.sets[read.sets == sample]
    if (length(input.reads) == 0) {
      print(paste0(sample, ": no read directory found. Skipping."))
      next
    }
    set.reads = list.files(paste0(actual.read.dir, "/", input.reads),
                           full.names = TRUE)

    ##########
    # Build HISAT2 index from missing targets (per-sample, unavoidable)
    index.dir = paste0(species.dir, "/index")
    if (!dir.exists(index.dir)) dir.create(index.dir)
    system(paste0(hisat2.path, "hisat2-build -f ", missing.ref, " ",
                  index.dir, "/reference"),
           ignore.stdout = quiet, ignore.stderr = quiet)

    # Map reads with permissive settings for divergent sequences
    system(paste0(hisat2.path, "hisat2 -q -x ", index.dir, "/reference",
                  " -1 ", set.reads[1], " -2 ", set.reads[2],
                  " -S ", species.dir, "/mapped_reads.sam",
                  " --mp 1,0 --sp 1,0 --score-min L,0.0,-0.3",
                  " --threads ", threads),
           ignore.stdout = quiet, ignore.stderr = quiet)

    # Extract: both mapped, R1-mapped/R2-unmapped, R2-mapped/R1-unmapped
    system(paste0(samtools.path, "samtools view -@ ", threads,
                  " -b -F 4 ", species.dir, "/mapped_reads.sam",
                  " > ", species.dir, "/mapped_all.bam"),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(samtools.path, "samtools view -@ ", threads,
                  " -b -f 4 -F 264 ", species.dir, "/mapped_reads.sam",
                  " > ", species.dir, "/mapped1.bam"),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(samtools.path, "samtools view -@ ", threads,
                  " -b -f 8 -F 260 ", species.dir, "/mapped_reads.sam",
                  " > ", species.dir, "/mapped2.bam"),
           ignore.stdout = quiet, ignore.stderr = quiet)

    system(paste0(samtools.path, "samtools merge -f -@ ", threads, " ",
                  species.dir, "/mapped_combined.bam ",
                  species.dir, "/mapped_all.bam ",
                  species.dir, "/mapped1.bam ",
                  species.dir, "/mapped2.bam"),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(samtools.path, "samtools sort -n -@ ", threads, " ",
                  species.dir, "/mapped_combined.bam",
                  " -o ", species.dir, "/mapped_sort.bam"),
           ignore.stdout = quiet, ignore.stderr = quiet)

    system(paste0("rm -f ",
                  species.dir, "/mapped_reads.sam ",
                  species.dir, "/mapped_all.bam ",
                  species.dir, "/mapped1.bam ",
                  species.dir, "/mapped2.bam ",
                  species.dir, "/mapped_combined.bam"))

    # Extract FASTQ from sorted BAM
    reads.dir = paste0(species.dir, "/temp_reads/sample")
    dir.create(reads.dir, recursive = TRUE, showWarnings = FALSE)
    system(paste0(samtools.path, "samtools fastq -@ ", threads, " ",
                  species.dir, "/mapped_sort.bam",
                  " -1 ", reads.dir, "/sample_READ1.fastq.gz",
                  " -2 ", reads.dir, "/sample_READ2.fastq.gz"),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0("rm -f ", species.dir, "/mapped_sort.bam"))

    # Merge paired-end reads with fastp before assembly
    PhyloProcessR::mergePairedEndReads(
      input.reads      = paste0(species.dir, "/temp_reads"),
      output.directory = paste0(species.dir, "/temp-merged-reads"),
      fastp.path       = fastp.path,
      threads          = threads,
      mem              = memory,
      overwrite        = TRUE,
      quiet            = quiet
    )

    merged.file = paste0(species.dir,
                         "/temp-merged-reads/sample/sample_READ3.fastq.gz")
    if (file.exists(merged.file) && file.size(merged.file) == 0) {
      system(paste0("rm -f ", merged.file))
    }
    system(paste0("rm -rf ", species.dir, "/temp_reads"))

    temp.read.path = list.files(paste0(species.dir, "/temp-merged-reads/sample"),
                                full.names = TRUE)

    # De novo assembly with SPAdes
    spades.contigs = PhyloProcessR::runSpades(
      read.paths          = temp.read.path,
      full.path.spades    = spades.path,
      mismatch.corrector  = TRUE,
      isolate             = FALSE,
      quiet               = quiet,
      read.contigs        = TRUE,
      threads             = threads,
      memory              = memory
    )

    system(paste0("rm -rf ",
                  species.dir, "/temp-merged-reads ",
                  index.dir))

    spades.contigs = spades.contigs[Biostrings::width(spades.contigs) >= 100]

    if (length(spades.contigs) == 0) {
      print(paste0(sample, ": SPAdes produced no contigs — saving original matches only."))
      found.contigs = Biostrings::readDNAStringSet(
        paste0(species.dir, "/", sample, "_matching-contigs.fa"))
      Biostrings::writeXStringSet(found.contigs, out.file)
      next
    }

    # BLAST new contigs against the per-sample missing reference
    contigs.file = paste0(species.dir, "/blast_contigs.fa")
    Biostrings::writeXStringSet(spades.contigs, contigs.file)

    system(paste0(blast.path, "makeblastdb -in ", missing.ref,
                  " -parse_seqids -dbtype nucl -out ", species.dir, "/blast_db"),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(blast.path, "blastn -task dc-megablast",
                  " -db ", species.dir, "/blast_db",
                  " -query ", contigs.file,
                  " -out ", species.dir, "/blast_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend",
                  " sstart send evalue bitscore qlen slen gaps\"",
                  " -num_threads ", threads),
           ignore.stdout = quiet, ignore.stderr = quiet)

    system(paste0("rm -f ",
                  species.dir, "/blast_db* ",
                  contigs.file))

    blast.out2 = paste0(species.dir, "/blast_match.txt")
    if (!file.exists(blast.out2) || file.size(blast.out2) == 0) {
      print(paste0(sample, ": no new contigs matched missing targets — saving originals."))
      found.contigs = Biostrings::readDNAStringSet(
        paste0(species.dir, "/", sample, "_matching-contigs.fa"))
      Biostrings::writeXStringSet(found.contigs, out.file)
      system(paste0("rm -f ", blast.out2))
      next
    }

    match.data2 = data.table::fread(blast.out2, sep = "\t", header = FALSE,
                                    stringsAsFactors = FALSE)
    system(paste0("rm -f ", blast.out2))
    data.table::setnames(match.data2, headers)

    filt.data2 = match.data2[match.data2$matches > min.match.length, ]
    filt.data2 = filt.data2[filt.data2$pident >= min.match.percent, ]
    filt.data2 = filt.data2[filt.data2$matches >= ((min.match.coverage / 100) * filt.data2$tLen), ]

    if (nrow(filt.data2) == 0) {
      print(paste0(sample, ": no new contigs passed filters — saving originals."))
      found.contigs = Biostrings::readDNAStringSet(
        paste0(species.dir, "/", sample, "_matching-contigs.fa"))
      Biostrings::writeXStringSet(found.contigs, out.file)
      next
    }

    data.table::setorder(filt.data2, qName, tName, -pident, -bitscore, evalue)
    new.contigs = spades.contigs[names(spades.contigs) %in% filt.data2$qName]

    write.table(filt.data2,
                file = paste0(species.dir, "/found-missing-blast-match.txt"),
                row.names = FALSE, quote = FALSE, sep = "\t")

    # Combine newly assembled + original matching contigs
    found.contigs  = Biostrings::readDNAStringSet(
      paste0(species.dir, "/", sample, "_matching-contigs.fa"))
    output.contigs = append(new.contigs, found.contigs)

    Biostrings::writeXStringSet(output.contigs, out.file)

    print(paste0(sample, " Phase 2 complete: ",
                 length(unique(found.data$tName)), " original targets + ",
                 length(new.contigs), " newly assembled."))

    rm(spades.contigs, match.data2, filt.data2, new.contigs,
       found.contigs, output.contigs)
    gc()

  }, error = function(e) {
    warning(file.names[i], " Phase 2 failed: ", conditionMessage(e))
  })
  } # end Phase 2 loop

  # Clean up shared reference BLAST DB
  system(paste0("rm -rf ", db.dir))

} # end function
