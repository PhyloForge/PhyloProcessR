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
#'   Must be paired (non-merged) reads. Common choices:
#'   \code{"decontaminated-reads"} (default), \code{"cleaned-reads"}. All read
#'   pairs in the sample folder are mapped, so multi-lane samples are handled.
#'   Merged reads cannot be used because HISAT2 expects paired input.
#'
#' @param reference path to the FASTA file of reference/target sequences used
#'   for BLAST matching and HISAT2 mapping.
#'
#' @param output.directory path to the directory where the per-sample working
#'   files are written. Default: \code{"expand-missing-assembly"}.
#'
#' @param expanded.directory path to the directory where the final expanded
#'   assemblies are written, one \code{<sample>.fa} file per sample. If
#'   \code{NULL} they are written to
#'   \code{output.directory/expanded-assemblies/}. Default: \code{NULL}.
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
#'   target from across all other samples -- a closer match to the true sequence,
#'   giving better read recovery. \code{"reference"} uses the original probe/
#'   bait sequence from the reference file -- useful when cross-sample contigs
#'   are unavailable or of low quality.
#'
#' @param recover.all.missing logical; if \code{TRUE}, a second recovery pass
#'   is performed for loci that are absent from every sample's assembly (i.e.
#'   not present in \code{final.save} at all). These loci are always mapped
#'   against the original reference sequences since no cross-sample contig
#'   exists. This pass can be time-consuming for large datasets. Default:
#'   \code{FALSE}.
#'
#' @param mismatch.corrector logical; if \code{TRUE} passes \code{--careful} to
#'   SPAdes when the recovered reads are assembled. Default: \code{TRUE}.
#'
#' @param kmer.values integer vector of k-mer sizes passed to SPAdes when the
#'   recovered reads are assembled. Default: \code{c(33, 55, 77, 99, 127)}.
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
#'   \code{expanded.directory/<sample>.fa}. Newly assembled
#'   contigs are renamed to their matched target, so the expanded assembly uses
#'   the same contig names as the input assembly.
#'
#' @export

expandMissingAssembly = function(assembly.directory = NULL,
                                 read.directory = NULL,
                                 mapping.reads = "decontaminated-reads",
                                 reference = NULL,
                                 output.directory = "expand-missing-assembly",
                                 expanded.directory = NULL,
                                 min.match.percent = 60,
                                 min.match.length = 100,
                                 min.match.coverage = 35,
                                 phase2.reference = c("contig", "reference"),
                                 recover.all.missing = FALSE,
                                 mismatch.corrector = TRUE,
                                 kmer.values = c(33, 55, 77, 99, 127),
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

  # Output directories. The finished assemblies can be placed outside the working
  # directory, so they can join the numbered contig directories.
  expanded.dir = if (is.null(expanded.directory)) {
    paste0(output.directory, "/expanded-assemblies")
  } else {
    expanded.directory
  }
  if (!dir.exists(output.directory)) dir.create(output.directory, recursive = TRUE)
  if (!dir.exists(expanded.dir))     dir.create(expanded.dir, recursive = TRUE)

  # Only FASTA files are samples. Hidden files and other output are ignored.
  fasta.pattern = "\\.fa$|\\.fas$|\\.fasta$|\\.fna$"
  file.names = list.files(assembly.directory, pattern = fasta.pattern)
  read.sets  = list.files(actual.read.dir)

  if (length(file.names) == 0) stop("No assembly files found in assembly.directory.")

  # Writes the contigs that already matched a target. Used whenever no new contig
  # is added for a sample.
  saveOriginals = function(sample.name, sample.dir, save.file) {
    found.contigs = Biostrings::readDNAStringSet(
      paste0(sample.dir, "/", sample.name, "_matching-contigs.fa"))
    Biostrings::writeXStringSet(found.contigs, save.file)
  }

  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  ###############################################################################
  ## Build BLAST reference DB once (major speedup -- was rebuilt per sample)
  ###############################################################################
  db.dir = paste0(output.directory, "/blast_ref_db")
  if (!dir.exists(db.dir)) dir.create(db.dir)
  system(paste0(blast.path, "makeblastdb -in ", shQuote(reference),
                " -parse_seqids -dbtype nucl -out ", shQuote(paste0(db.dir, "/nucl-blast_db"))),
         ignore.stdout = quiet, ignore.stderr = quiet)

  ###############################################################################
  ## Phase 1: BLAST each assembly against reference -- parallelized
  ###############################################################################
  parallel::mclapply(seq_along(file.names), function(i) {
  tryCatch({

    sample      = gsub(fasta.pattern, "", file.names[i])
    species.dir = paste0(output.directory, "/", sample)
    if (!dir.exists(species.dir)) dir.create(species.dir, recursive = TRUE)

    # Skip if already done
    if (overwrite == FALSE &&
        file.exists(paste0(species.dir, "/filtered-blast-match.txt"))) {
      print(paste0(sample, " Phase 1 already done, skipping."))
      return(NULL)
    }

    # BLAST assembly against shared reference DB
    blast.out = paste0(species.dir, "/target-blast-match.txt")
    system(paste0(
      blast.path, "blastn -task dc-megablast -db ", shQuote(paste0(db.dir, "/nucl-blast_db")),
      " -evalue 0.001",
      " -query ", shQuote(paste0(assembly.directory, "/", file.names[i])),
      " -out ", shQuote(blast.out),
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
    unlink(blast.out)

    # Filter by identity, length, and coverage
    filt.data = match.data[match.data$matches > min.match.length, ]
    filt.data = filt.data[filt.data$pident >= min.match.percent, ]
    filt.data = filt.data[filt.data$matches >= ((min.match.coverage / 100) * filt.data$tLen), ]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, ": no matches passed filters. Skipping."))
      return(NULL)
    }

    # Bitscore comes before identity so a long strong hit beats a short exact one
    data.table::setorder(filt.data, qName, tName, -bitscore, -pident, evalue)

    sample.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", file.names[i]))
    # BLAST reports the first word of the header, so the names are cut to match
    names(sample.contigs) = gsub(" .*", "", names(sample.contigs))
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
  # Each sample is collected once and joined at the end. Growing the set one
  # contig at a time is quadratic on large datasets.
  match.list = vector("list", length(file.names))

  for (i in seq_along(file.names)) {
    sample      = gsub(fasta.pattern, "", file.names[i])
    species.dir = paste0(output.directory, "/", sample)
    match.file  = paste0(species.dir, "/", sample, "_matching-contigs.fa")
    blast.file  = paste0(species.dir, "/filtered-blast-match.txt")

    if (!file.exists(match.file) || !file.exists(blast.file)) next

    match.contigs = Biostrings::readDNAStringSet(match.file)
    filt.data     = read.table(blast.file, sep = "\t", header = TRUE,
                               stringsAsFactors = FALSE)

    keep.index = match(filt.data$qName, names(match.contigs))
    filt.data  = filt.data[!is.na(keep.index), ]
    keep.index = keep.index[!is.na(keep.index)]
    if (length(keep.index) == 0) next

    temp.match = match.contigs[keep.index]
    names(temp.match) = filt.data$tName
    match.list[[i]] = temp.match
  }

  match.list  = match.list[!vapply(match.list, is.null, logical(1))]
  all.matches = if (length(match.list) == 0) Biostrings::DNAStringSet() else do.call(c, match.list)

  # Keep longest contig per target name
  if (length(all.matches) > 0) {
    all.matches = all.matches[order(names(all.matches), -Biostrings::width(all.matches))]
    final.save  = all.matches[!duplicated(names(all.matches))]
  } else {
    final.save = all.matches
  }

  Biostrings::writeXStringSet(final.save,
                              paste0(output.directory, "/unique_matches.fa"))

  if (length(final.save) == 0 && recover.all.missing == FALSE) {
    message("No cross-sample matches found. Cannot proceed to Phase 2.")
    unlink(db.dir, recursive = TRUE)
    return(invisible(NULL))
  }

  # Load reference sequences if needed for phase2.reference="reference" or
  # recover.all.missing=TRUE (lazy -- only pay the I/O cost when required)
  if (phase2.reference == "reference" || recover.all.missing == TRUE) {
    reference.seqs = Biostrings::readDNAStringSet(reference)
  }

  ###############################################################################
  ## Phase 2: Map reads to missing targets, assemble, expand -- sequential so
  ##           each tool (HISAT2, SPAdes, samtools) gets the full thread budget
  ###############################################################################
  for (i in seq_along(file.names)) {
  tryCatch({

    sample      = gsub(fasta.pattern, "", file.names[i])
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
      print(paste0(sample, ": all targets already assembled -- saving originals."))
      saveOriginals(sample, species.dir, out.file)
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
    set.reads = set.reads[grep("fastq|fq", basename(set.reads))]

    # Reads are selected by name, not by position, so every lane is mapped and an
    # extra file in the folder cannot shift the pair.
    read1 = sort(set.reads[grep("_1\\.f|-1\\.f|_R1[_.-]|-R1[_.-]|READ1", basename(set.reads))])
    read2 = sort(set.reads[grep("_2\\.f|-2\\.f|_R2[_.-]|-R2[_.-]|READ2", basename(set.reads))])

    if (length(read1) == 0 || length(read1) != length(read2)) {
      warning(sample, ": found ", length(read1), " read1 and ", length(read2),
              " read2 files in ", mapping.reads, ". Paired reads are required. Skipping.")
      next
    }

    ##########
    # Build HISAT2 index from missing targets (per-sample, unavoidable)
    index.dir = paste0(species.dir, "/index")
    if (!dir.exists(index.dir)) dir.create(index.dir, recursive = TRUE)
    system(paste0(hisat2.path, "hisat2-build -f ", shQuote(missing.ref), " ",
                  shQuote(paste0(index.dir, "/reference"))),
           ignore.stdout = quiet, ignore.stderr = quiet)

    # Map reads with permissive settings for divergent sequences. HISAT2 takes a
    # comma separated list, so all lanes are mapped in one run.
    system(paste0(hisat2.path, "hisat2 -q -x ", shQuote(paste0(index.dir, "/reference")),
                  " -1 ", paste(shQuote(read1), collapse = ","),
                  " -2 ", paste(shQuote(read2), collapse = ","),
                  " -S ", shQuote(paste0(species.dir, "/mapped_reads.sam")),
                  " --mp 1,0 --sp 1,0 --score-min L,0.0,-0.3",
                  " --threads ", threads),
           ignore.stdout = quiet, ignore.stderr = quiet)

    # Extract: both mapped, R1-mapped/R2-unmapped, R2-mapped/R1-unmapped
    sam.file = paste0(species.dir, "/mapped_reads.sam")
    system(paste0(samtools.path, "samtools view -@ ", threads,
                  " -b -F 4 ", shQuote(sam.file),
                  " > ", shQuote(paste0(species.dir, "/mapped_all.bam"))),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(samtools.path, "samtools view -@ ", threads,
                  " -b -f 4 -F 264 ", shQuote(sam.file),
                  " > ", shQuote(paste0(species.dir, "/mapped1.bam"))),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(samtools.path, "samtools view -@ ", threads,
                  " -b -f 8 -F 260 ", shQuote(sam.file),
                  " > ", shQuote(paste0(species.dir, "/mapped2.bam"))),
           ignore.stdout = quiet, ignore.stderr = quiet)

    system(paste0(samtools.path, "samtools merge -f -@ ", threads, " ",
                  shQuote(paste0(species.dir, "/mapped_combined.bam")), " ",
                  shQuote(paste0(species.dir, "/mapped_all.bam")), " ",
                  shQuote(paste0(species.dir, "/mapped1.bam")), " ",
                  shQuote(paste0(species.dir, "/mapped2.bam"))),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(samtools.path, "samtools sort -n -@ ", threads, " ",
                  shQuote(paste0(species.dir, "/mapped_combined.bam")),
                  " -o ", shQuote(paste0(species.dir, "/mapped_sort.bam"))),
           ignore.stdout = quiet, ignore.stderr = quiet)

    unlink(paste0(species.dir, c("/mapped_reads.sam", "/mapped_all.bam",
                                 "/mapped1.bam", "/mapped2.bam",
                                 "/mapped_combined.bam")))

    # Extract FASTQ from sorted BAM
    reads.dir = paste0(species.dir, "/temp_reads/sample")
    dir.create(reads.dir, recursive = TRUE, showWarnings = FALSE)
    system(paste0(samtools.path, "samtools fastq -@ ", threads, " ",
                  shQuote(paste0(species.dir, "/mapped_sort.bam")),
                  " -1 ", shQuote(paste0(reads.dir, "/sample_READ1.fastq.gz")),
                  " -2 ", shQuote(paste0(reads.dir, "/sample_READ2.fastq.gz"))),
           ignore.stdout = quiet, ignore.stderr = quiet)
    unlink(paste0(species.dir, "/mapped_sort.bam"))

    # An empty gzip file is about 20 bytes. Stop here when no read mapped, rather
    # than passing empty files to fastp and SPAdes.
    mapped.reads = paste0(reads.dir, c("/sample_READ1.fastq.gz", "/sample_READ2.fastq.gz"))
    if (any(file.exists(mapped.reads) == FALSE) || min(file.size(mapped.reads)) < 100) {
      print(paste0(sample, ": no reads mapped to the missing targets -- saving originals."))
      saveOriginals(sample, species.dir, out.file)
      unlink(paste0(species.dir, "/temp_reads"), recursive = TRUE)
      unlink(index.dir, recursive = TRUE)
      next
    }

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
      unlink(merged.file)
    }
    unlink(paste0(species.dir, "/temp_reads"), recursive = TRUE)

    temp.read.path = list.files(paste0(species.dir, "/temp-merged-reads/sample"),
                                full.names = TRUE)

    # De novo assembly with SPAdes
    spades.contigs = PhyloProcessR::runSpades(
      read.paths          = temp.read.path,
      full.path.spades    = spades.path,
      mismatch.corrector  = mismatch.corrector,
      isolate             = FALSE,
      kmer.values         = kmer.values,
      quiet               = quiet,
      read.contigs        = TRUE,
      clean               = TRUE,
      threads             = threads,
      memory              = memory
    )

    unlink(c(paste0(species.dir, "/temp-merged-reads"), index.dir), recursive = TRUE)

    spades.contigs = spades.contigs[Biostrings::width(spades.contigs) >= 100]

    if (length(spades.contigs) == 0) {
      print(paste0(sample, ": SPAdes produced no contigs -- saving original matches only."))
      saveOriginals(sample, species.dir, out.file)
      next
    }

    # BLAST new contigs against the per-sample missing reference
    contigs.file = paste0(species.dir, "/blast_contigs.fa")
    Biostrings::writeXStringSet(spades.contigs, contigs.file)

    system(paste0(blast.path, "makeblastdb -in ", shQuote(missing.ref),
                  " -parse_seqids -dbtype nucl -out ", shQuote(paste0(species.dir, "/blast_db"))),
           ignore.stdout = quiet, ignore.stderr = quiet)
    system(paste0(blast.path, "blastn -task dc-megablast",
                  " -db ", shQuote(paste0(species.dir, "/blast_db")),
                  " -query ", shQuote(contigs.file),
                  " -out ", shQuote(paste0(species.dir, "/blast_match.txt")),
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend",
                  " sstart send evalue bitscore qlen slen gaps\"",
                  " -num_threads ", threads),
           ignore.stdout = quiet, ignore.stderr = quiet)

    unlink(c(list.files(species.dir, pattern = "^blast_db", full.names = TRUE), contigs.file))

    blast.out2 = paste0(species.dir, "/blast_match.txt")
    if (!file.exists(blast.out2) || file.size(blast.out2) == 0) {
      print(paste0(sample, ": no new contigs matched missing targets -- saving originals."))
      saveOriginals(sample, species.dir, out.file)
      unlink(blast.out2)
      next
    }

    match.data2 = data.table::fread(blast.out2, sep = "\t", header = FALSE,
                                    stringsAsFactors = FALSE)
    unlink(blast.out2)
    data.table::setnames(match.data2, headers)

    filt.data2 = match.data2[match.data2$matches > min.match.length, ]
    filt.data2 = filt.data2[filt.data2$pident >= min.match.percent, ]
    filt.data2 = filt.data2[filt.data2$matches >= ((min.match.coverage / 100) * filt.data2$tLen), ]

    if (nrow(filt.data2) == 0) {
      print(paste0(sample, ": no new contigs passed filters -- saving originals."))
      saveOriginals(sample, species.dir, out.file)
      next
    }

    # Bitscore comes before identity so a long strong hit beats a short exact one
    data.table::setorder(filt.data2, qName, tName, -bitscore, -pident, evalue)

    # Keeps the best hit for each new contig and renames the contig to its target.
    # The expanded assembly then uses the same contig names as the input assembly.
    best.data2 = filt.data2[duplicated(filt.data2$qName) == FALSE, ]
    keep.index = match(best.data2$qName, names(spades.contigs))
    best.data2 = best.data2[!is.na(keep.index), ]
    keep.index = keep.index[!is.na(keep.index)]

    if (length(keep.index) == 0) {
      print(paste0(sample, ": no new contigs passed filters -- saving originals."))
      saveOriginals(sample, species.dir, out.file)
      next
    }

    new.contigs = spades.contigs[keep.index]
    names(new.contigs) = best.data2$tName

    write.table(filt.data2,
                file = paste0(species.dir, "/found-missing-blast-match.txt"),
                row.names = FALSE, quote = FALSE, sep = "\t")

    # Combine original matching contigs + newly assembled. The original contigs
    # come first so they keep their name if a target is recovered twice.
    found.contigs  = Biostrings::readDNAStringSet(
      paste0(species.dir, "/", sample, "_matching-contigs.fa"))
    output.contigs = append(found.contigs, new.contigs)
    names(output.contigs) = make.unique(names(output.contigs), sep = "_")

    Biostrings::writeXStringSet(output.contigs, out.file)

    print(paste0(sample, " Phase 2 complete: ",
                 length(unique(found.data$tName)), " original targets + ",
                 length(new.contigs), " newly assembled."))

    rm(spades.contigs, match.data2, filt.data2, best.data2, new.contigs,
       found.contigs, output.contigs)
    gc()

  }, error = function(e) {
    warning(file.names[i], " Phase 2 failed: ", conditionMessage(e))
  })
  } # end Phase 2 loop

  # Clean up shared reference BLAST DB
  unlink(db.dir, recursive = TRUE)

} # end function
