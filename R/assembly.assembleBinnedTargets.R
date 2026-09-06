#' @title assembleBinnedTargets
#'
#' @description Assembles each target locus on its own, to recover targets that
#'   the whole-library assembly lost and to extend the targets it found. The
#'   function maps the paired reads of a sample to a set of bait sequences, puts
#'   each read pair in the bin of its best-scoring target, and then runs a
#'   separate megahit assembly for every bin. A pair enters a bin when either
#'   mate aligns, so the mate of an anchored read adds sequence beyond the
#'   target. Each locus also gets its own coverage distribution, which protects a
#'   low-coverage locus from the coverage cutoffs of a whole-library assembly.
#'
#'   The step is additive. It does not replace \code{assembleSpades}. Give it the
#'   assembly directory of the main pipeline and it keeps the longer sequence for
#'   each target, so a locus that assembled well keeps its own contig.
#'
#'   Divergence is handled at the identification step, not at the mapping step. A
#'   target locus assembles in the draft assembly whatever its divergence,
#'   because de novo assembly needs no reference. The locus is then lost when
#'   \code{blastn} cannot match the contig to a probe. Give
#'   \code{draft.assembly.directory} and LAST searches the draft assembly of the
#'   sample against the targets. LAST trains a scoring matrix on the sample, so
#'   it finds contigs at 35 percent divergence that \code{blastn} misses. Those
#'   contigs then become the baits, and bwa maps the reads to sequence of their
#'   own sample at full identity.
#'
#' @details Set \code{locus.set = "missing"} to work only on the targets a sample
#'   does not have. This is much faster and it is the safe first run. Set
#'   \code{locus.set = "all"} to also extend the targets the sample has.
#'
#'   The per-bin assemblies control the run time. megahit takes about 1.1 seconds
#'   for a bin on one thread, at any bin depth. A sample with 14,000 bins
#'   therefore costs about 4.5 CPU hours, or about 20 minutes on 16 threads.
#'   Timing numbers are in HANDOFF.md.
#'
#'   The function deletes each bin directory as soon as it reads the contigs. A
#'   megahit run on a small bin writes about 27 files. Tens of thousands of bins
#'   can therefore use all the inodes of a shared filesystem.
#'
#' @param read.directory path to the top-level processed reads directory (for
#'   example \code{"processed-reads"}).
#'
#' @param mapping.reads name of the subdirectory in \code{read.directory} that
#'   holds the per-sample read folders. The reads must be paired. Merged reads
#'   cannot be used, because the mate of an anchored read supplies the flanking
#'   sequence. Default: \code{"decontaminated-reads"}.
#'
#' @param target.markers path to the FASTA file of target markers. Each sequence
#'   becomes one bin.
#'
#' @param assembly.directory path to the per-sample assemblies of the main
#'   pipeline, one \code{<sample>.fa} file per sample, with the contigs named
#'   after their target. The new contigs are merged with these. If \code{NULL}
#'   the output holds only the new contigs and \code{locus.set} must be
#'   \code{"all"}. Default: \code{NULL}.
#'
#' @param draft.assembly.directory path to the draft assemblies of the sample,
#'   one \code{<sample>.fa} file per sample, before the target filter removed
#'   the unmatched contigs (for example \code{2_reduced-redundancy}). LAST
#'   searches these contigs against the targets to find the divergent loci that
#'   \code{blastn} could not assign. If \code{NULL} this search is skipped and
#'   LAST is not needed. Default: \code{NULL}.
#'
#' @param output.directory path to the working directory. Default:
#'   \code{"binned-target-assembly"}.
#'
#' @param binned.directory path to the directory for the final per-sample FASTA
#'   files. If \code{NULL} they go to \code{output.directory/binned-assemblies}.
#'   Default: \code{NULL}.
#'
#' @param locus.set which targets to bin. \code{"missing"} uses only the targets
#'   absent from the assembly of that sample. \code{"all"} uses every target, so
#'   the targets the sample already has are extended as well. Default:
#'   \code{"all"}.
#'
#' @param bait.source what to map the reads against in the first round.
#'   \code{"hybrid"} uses a sequence of the sample when there is one that covers
#'   \code{min.bait.coverage} of the target, and the reference for the other
#'   targets. Reads map to sequence of their own sample at full identity, so bwa
#'   stays sensitive. \code{"reference"} always uses the target markers.
#'   Default: \code{"hybrid"}.
#'
#' @param min.bait.coverage least part of the target length that a sequence of
#'   the sample must cover before it is used as the bait, on a scale of 0 to 1.
#'   A short fragment recruits reads only across itself, so the reference is used
#'   instead. Default: \code{0.5}.
#'
#' @param rescue.missing logical. \code{TRUE} recruits reads with LAST for the
#'   targets that have no sequence in this sample, assembles them, and uses the
#'   result as the bait. bwa cannot recruit a read that is 35 percent divergent
#'   from its bait, so without this step a divergent target that the draft
#'   assembly also lost cannot be recovered. Default: \code{TRUE}.
#'
#' @param iterations number of bait-and-assemble rounds. One round recovers about
#'   one insert length of flank on each side. A later round baits with the
#'   contigs of the round before it and adds about one more insert length, at a
#'   growing risk of extension into a repeat. Rounds after the first only process
#'   the targets that produced a contig. Default: \code{1}.
#'
#' @param min.pairs minimum read pairs a bin must hold before it is assembled.
#'   Default: \code{6}.
#'
#' @param max.pairs maximum read pairs kept per bin. A very deep locus costs
#'   assembly time and gives the assembler no more information. \code{0} removes the
#'   limit. Default: \code{3000}.
#'
#' @param min.contig.length minimum length in base pairs of an assembled contig.
#'   Default: \code{100}.
#'
#' @param min.match.percent minimum BLAST percent identity of a new contig
#'   against its target. Default: \code{60}.
#'
#' @param min.match.length minimum BLAST alignment length in base pairs.
#'   Default: \code{60}.
#'
#' @param min.match.coverage minimum percentage of the target length that the
#'   BLAST match must cover. This is lower than the whole-library default,
#'   because a binned contig is anchored to its own target already. Default:
#'   \code{20}.
#'
#' @param max.extension maximum base pairs a contig can add to each side of its
#'   bait in one round. This stops a round that extends into a transposable
#'   element. The length rule would otherwise keep that contig. Default:
#'   \code{1000}.
#'
#' @param max.target.hits maximum number of targets a new contig may match. A
#'   contig above this is a shared repeat and is dropped. Default: \code{5}.
#'
#' @param multi.copy what to do with a target that has more than one contig in
#'   \code{assembly.directory}. \code{"keep"} leaves those targets unchanged. One
#'   binned contig cannot represent two copies. \code{"longest"} treats them like
#'   any other target and keeps the longest sequence. \code{reduceRedundancy}
#'   collapses the near-identical haplotype contigs before this step, so most
#'   targets with two contigs here are paralogs. Default: \code{"keep"}.
#'
#' @param kmer.values integer vector of k-mer sizes for the assembly. Every value
#'   must be below the read length. megahit also needs odd values between 15 and
#'   255, and a step of 28 or less between them. Fewer values are faster.
#'   Default: \code{c(21, 33, 55, 77, 99)}.
#'
#' @param memory RAM in GB for the step. Each parallel assembly job gets a share.
#'   Default: \code{8}.
#'
#' @param threads number of CPU threads. The mapping step uses all of them. The
#'   assembly step runs this many single-threaded megahit jobs at once. Default:
#'   \code{1}.
#'
#' @param bwa.path path to the directory that holds \code{bwa}. If \code{NULL}
#'   the program must be on the system PATH. Default: \code{NULL}.
#'
#' @param samtools.path path to the directory that holds \code{samtools}.
#'   Default: \code{NULL}.
#'
#' @param megahit.path path to the directory that holds \code{megahit}.
#'   Default: \code{NULL}.
#'
#' @param spades.path path to the directory that holds \code{spades.py}. Only
#'   \code{rescue.missing = TRUE} needs it. Default: \code{NULL}.
#'
#' @param last.path path to the directory that holds \code{lastdb},
#'   \code{lastal} and \code{maf-convert}. Default: \code{NULL}.
#'
#' @param overwrite logical. \code{TRUE} runs a sample again when its output
#'   exists. Default: \code{FALSE}.
#'
#' @param quiet logical. \code{TRUE} hides the output of the external programs.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. The function writes one \code{<sample>.fa}
#'   file per sample to \code{binned.directory}, and one
#'   \code{<sample>_binned-stats.txt} table to the working directory of the
#'   sample. The table gives the bait source, the reads in the bin, the new
#'   length and the previous length for every target.
#'
#' @export

assembleBinnedTargets = function(read.directory = NULL,
                                 mapping.reads = "decontaminated-reads",
                                 target.markers = NULL,
                                 assembly.directory = NULL,
                                 draft.assembly.directory = NULL,
                                 output.directory = "binned-target-assembly",
                                 binned.directory = NULL,
                                 locus.set = c("all", "missing"),
                                 bait.source = c("hybrid", "reference"),
                                 min.bait.coverage = 0.5,
                                 rescue.missing = TRUE,
                                 iterations = 1,
                                 min.pairs = 6,
                                 max.pairs = 3000,
                                 min.contig.length = 100,
                                 min.match.percent = 60,
                                 min.match.length = 60,
                                 min.match.coverage = 20,
                                 max.extension = 1000,
                                 max.target.hits = 5,
                                 multi.copy = c("keep", "longest"),
                                 kmer.values = c(21, 33, 55, 77, 99),
                                 memory = 8,
                                 threads = 1,
                                 bwa.path = NULL,
                                 samtools.path = NULL,
                                 megahit.path = NULL,
                                 spades.path = NULL,
                                 last.path = NULL,
                                 overwrite = FALSE,
                                 quiet = TRUE) {

  locus.set   = match.arg(locus.set)
  bait.source = match.arg(bait.source)
  multi.copy  = match.arg(multi.copy)

  if (is.null(read.directory) == TRUE) {
    stop("A read directory is needed. Set read.directory.")
  }
  if (is.null(target.markers) == TRUE || file.exists(target.markers) == FALSE) {
    stop("The target marker file was not found. Set target.markers.")
  }
  if (locus.set == "missing" && is.null(assembly.directory) == TRUE) {
    stop("locus.set = \"missing\" needs an assembly.directory to know what is missing.")
  }
  if (iterations < 1) stop("iterations must be 1 or greater.")

  # megahit rejects an even k, a k outside 15 to 255, and a step above 28
  if (any(kmer.values %% 2 == 0) || any(kmer.values < 15) ||
      any(kmer.values > 255) || any(diff(sort(kmer.values)) > 28)) {
    stop("kmer.values must be odd, between 15 and 255, and no more than 28 apart.")
  }

  bwa.command      = .toolCommand("bwa", bwa.path)
  samtools.command = .toolCommand("samtools", samtools.path)
  megahit.command  = .toolCommand("megahit", megahit.path)

  # Only the rescue step still uses SPAdes, so it is the only thing that needs it
  use.rescue = rescue.missing == TRUE && bait.source != "reference"
  spades.command = NULL
  if (use.rescue == TRUE) spades.command = .toolCommand("spades.py", spades.path)

  # LAST does every divergent search in this function, so it is always needed
  use.draft = is.null(draft.assembly.directory) == FALSE
  lastdb.command     = .toolCommand("lastdb", last.path)
  lastal.command     = .toolCommand("lastal", last.path)
  mafconvert.command = .toolCommand("maf-convert", last.path)

  actual.read.dir = paste0(sub("/+$", "", read.directory), "/", mapping.reads)
  if (dir.exists(actual.read.dir) == FALSE) {
    stop("The read directory ", actual.read.dir, " was not found.")
  }

  output.directory = sub("/+$", "", output.directory)
  if (is.null(binned.directory) == TRUE) {
    binned.directory = paste0(output.directory, "/binned-assemblies")
  }
  binned.directory = sub("/+$", "", binned.directory)
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(binned.directory, recursive = TRUE, showWarnings = FALSE)

  reference.seqs = Biostrings::readDNAStringSet(target.markers)
  names(reference.seqs) = gsub(" .*", "", names(reference.seqs))
  reference.seqs = reference.seqs[duplicated(names(reference.seqs)) == FALSE]
  if (length(reference.seqs) == 0) stop("No sequences were read from target.markers.")
  locus.names = names(reference.seqs)

  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore",
              "qLen", "tLen", "gaps")

  # One LAST database of the targets, shared by every sample
  db.dir = paste0(output.directory, "/target-databases")
  dir.create(db.dir, showWarnings = FALSE)
  target.copy = paste0(db.dir, "/targets.fa")
  Biostrings::writeXStringSet(reference.seqs, target.copy)

  last.db = paste0(db.dir, "/last-targets")
  .lastBuildDB(reference.file = target.copy,
               db.prefix = last.db,
               lastdb.command = lastdb.command,
               threads = threads,
               quiet = quiet)

  sample.names = list.dirs(actual.read.dir, full.names = FALSE, recursive = FALSE)
  sample.names = sample.names[nchar(sample.names) > 0]
  if (length(sample.names) == 0) {
    stop("No sample folders were found in ", actual.read.dir, ".")
  }

  assembly.memory = max(2, floor(memory / max(1, threads)))

  for (i in seq_along(sample.names)) {
  tryCatch({

    sample     = sample.names[i]
    sample.dir = paste0(output.directory, "/", sample)
    out.file   = paste0(binned.directory, "/", sample, ".fa")

    if (overwrite == FALSE && file.exists(out.file) == TRUE) {
      print(paste0(sample, " already finished, skipping."))
      next
    }
    if (overwrite == TRUE) unlink(sample.dir, recursive = TRUE)
    dir.create(sample.dir, recursive = TRUE, showWarnings = FALSE)

    read.pair = .pairSampleReads(paste0(actual.read.dir, "/", sample))
    if (is.null(read.pair) == TRUE) {
      warning(sample, ": paired reads were not found in ", mapping.reads, ". Skipping.")
      next
    }

    # bwa mem takes one file per mate, so several lanes are joined first
    read1 = read.pair$read1
    read2 = read.pair$read2
    if (length(read1) > 1) {
      read1 = paste0(sample.dir, "/all_R1.fastq")
      read2 = paste0(sample.dir, "/all_R2.fastq")
      .runCommand(paste0("gzip -cdf ", paste(shQuote(read.pair$read1), collapse = " "),
                         " > ", shQuote(read1)),
                  quiet = quiet, task = "read joining", keep.stdout = TRUE)
      .runCommand(paste0("gzip -cdf ", paste(shQuote(read.pair$read2), collapse = " "),
                         " > ", shQuote(read2)),
                  quiet = quiet, task = "read joining", keep.stdout = TRUE)
    }

    old.contigs = NULL
    if (is.null(assembly.directory) == FALSE) {
      old.file = paste0(sub("/+$", "", assembly.directory), "/", sample, ".fa")
      if (file.exists(old.file) == TRUE) {
        old.contigs = Biostrings::readDNAStringSet(old.file)
        names(old.contigs) = gsub(" .*", "", names(old.contigs))
      } else {
        print(paste0(sample, ": no previous assembly found, baiting with the reference."))
      }
    }
    own.contigs = .longestPerTarget(old.contigs, locus.names)

    # LAST finds the divergent contigs that the nucleotide BLAST of the main
    # pipeline could not assign to a target
    draft.contigs = NULL
    if (use.draft == TRUE) {
      draft.file = paste0(sub("/+$", "", draft.assembly.directory), "/", sample, ".fa")
      if (file.exists(draft.file) == TRUE) {
        draft.contigs = .identifyDraftContigs(
          draft.file = draft.file,
          last.db = last.db,
          work.dir = sample.dir,
          lastal.command = lastal.command,
          headers = headers,
          min.match.percent = min.match.percent,
          min.match.length = min.match.length,
          min.match.coverage = min.match.coverage,
          threads = threads,
          quiet = quiet)
        print(paste0(sample, ": LAST assigned ",
                     if (is.null(draft.contigs)) 0 else length(draft.contigs),
                     " draft contigs to a target."))
      } else {
        print(paste0(sample, ": no draft assembly found at ", draft.file, "."))
      }
    }

    target.names = locus.names
    if (locus.set == "missing" && is.null(own.contigs) == FALSE) {
      target.names = locus.names[locus.names %in% names(own.contigs) == FALSE]
    }
    if (length(target.names) == 0) {
      print(paste0(sample, ": no target needs binning. The previous assembly is saved unchanged."))
      if (is.null(old.contigs) == FALSE) Biostrings::writeXStringSet(old.contigs, out.file)
      next
    }

    # Targets with no sequence in this sample at all. bwa cannot recruit for
    # these when the sample is divergent, so LAST recruits the reads instead.
    seed.contigs = NULL
    if (use.rescue == TRUE) {
      have.names = c(names(own.contigs), names(draft.contigs))
      rescue.names = target.names[target.names %in% have.names == FALSE]

      if (length(rescue.names) > 0) {
        seed.contigs = .rescueMissingTargets(
          missing.seqs = reference.seqs[rescue.names],
          read.files = c(read1, read2),
          work.dir = sample.dir,
          lastdb.command = lastdb.command,
          lastal.command = lastal.command,
          mafconvert.command = mafconvert.command,
          samtools.command = samtools.command,
          spades.command = spades.command,
          headers = headers,
          kmer.values = kmer.values,
          min.contig.length = min.contig.length,
          min.match.percent = min.match.percent,
          min.match.length = min.match.length,
          min.match.coverage = min.match.coverage,
          memory = memory,
          threads = threads,
          quiet = quiet)

        print(paste0(sample, ": LAST rescued ",
                     if (is.null(seed.contigs)) 0 else length(seed.contigs),
                     " of ", length(rescue.names), " targets with no sequence."))
      }
    }

    bait.set = .buildBaitTable(reference.seqs = reference.seqs,
                               own.contigs = own.contigs,
                               draft.contigs = draft.contigs,
                               seed.contigs = seed.contigs,
                               target.names = target.names,
                               bait.source = bait.source,
                               min.bait.coverage = min.bait.coverage)

    print(paste0(sample, ": ", nrow(bait.set$table), " baits. Sources: ",
                 paste(names(table(bait.set$table$source)),
                       table(bait.set$table$source), sep = "=", collapse = ", ")))

    best.contigs = NULL
    stats.table  = bait.set$table
    stats.table$reads = 0

    for (round in seq_len(iterations)) {

      round.dir = paste0(sample.dir, "/round", round)
      unlink(round.dir, recursive = TRUE)
      dir.create(round.dir, recursive = TRUE, showWarnings = FALSE)

      bait.file = paste0(round.dir, "/bait.fa")
      Biostrings::writeXStringSet(bait.set$seqs, bait.file)
      .runCommand(paste0(bwa.command, " index ", shQuote(bait.file)),
                  quiet = quiet, task = "bait indexing")

      # Map once, then sort and index. samtools reads each bin from the index
      # and writes the FASTQ, so no alignment record is parsed here. The -k 15
      # and -T 25 options relax bwa, which matters only for a reference bait.
      bam.file = paste0(round.dir, "/mapped.bam")
      .runCommand(paste0(bwa.command, " mem -t ", threads,
                         " -k 15 -B 3 -O 5 -T 25 ",
                         shQuote(bait.file), " ", shQuote(read1), " ", shQuote(read2),
                         " | ", samtools.command, " sort -@ ", threads,
                         " -T ", shQuote(paste0(round.dir, "/sorting")),
                         " -o ", shQuote(bam.file)),
                  quiet = quiet, task = "read mapping")
      .runCommand(paste0(samtools.command, " index ", shQuote(bam.file)),
                  quiet = quiet, task = "BAM indexing")

      unlink(paste0(bait.file, c(".amb", ".ann", ".bwt", ".pac", ".sa")))

      # idxstats gives the reads per bait in one pass, so an empty bin never
      # reaches the assembler. Column 4 holds the unmapped mates placed at that bait.
      idx.file = paste0(round.dir, "/idxstats.txt")
      .runCommand(paste0(samtools.command, " idxstats ", shQuote(bam.file),
                         " > ", shQuote(idx.file)),
                  quiet = quiet, task = "bin counting", keep.stdout = TRUE)

      idx.data = utils::read.table(idx.file, sep = "\t", header = FALSE,
                                   stringsAsFactors = FALSE)
      colnames(idx.data) = c("bait", "length", "mapped", "unmapped")
      idx.data = idx.data[idx.data$bait != "*", ]

      bait.reads = idx.data$mapped + idx.data$unmapped
      names(bait.reads) = idx.data$bait
      round.reads = as.numeric(bait.reads[bait.set$table$bait])
      round.reads[is.na(round.reads)] = 0
      if (round == 1) stats.table$reads = round.reads[match(stats.table$locus, bait.set$table$locus)]

      run.index = which(round.reads >= (min.pairs * 2))
      print(paste0(sample, " round ", round, ": ", length(run.index), " of ",
                   nrow(bait.set$table), " bins hold at least ", min.pairs,
                   " read pairs. Assembling."))

      if (length(run.index) == 0) break

      # One megahit job per bin, in parallel. A bin has a few hundred reads. Many
      # single-threaded jobs are therefore faster than one threaded job.
      bin.dir = paste0(round.dir, "/bins")
      dir.create(bin.dir, showWarnings = FALSE)

      contig.list = parallel::mclapply(run.index, function(j) {
        # samtools -s takes the seed in the integer part and the fraction to
        # keep after the point. A very deep bin costs time and adds nothing.
        subsample = NULL
        if (max.pairs > 0 && round.reads[j] > (max.pairs * 2)) {
          subsample = format(42 + ((max.pairs * 2) / round.reads[j]), nsmall = 4)
        }
        .assembleOneBin(bait = bait.set$table$bait[j],
                        bam.file = bam.file,
                        bin.dir = bin.dir,
                        samtools.command = samtools.command,
                        megahit.command = megahit.command,
                        kmer.values = kmer.values,
                        subsample = subsample,
                        memory = assembly.memory,
                        min.contig.length = min.contig.length,
                        quiet = quiet)
      }, mc.cores = min(threads, length(run.index)), mc.preschedule = FALSE)

      contig.list = contig.list[vapply(contig.list, length, integer(1)) > 0]
      unlink(c(bin.dir, bam.file, paste0(bam.file, ".bai")), recursive = TRUE)

      if (length(contig.list) == 0) {
        print(paste0(sample, " round ", round, ": the assembly produced no contig."))
        break
      }

      round.best = .filterBinnedContigs(
        contigs = do.call(c, unname(contig.list)),
        bait.table = bait.set$table,
        bait.widths = stats::setNames(Biostrings::width(bait.set$seqs),
                                      bait.set$table$locus),
        db.prefix = last.db,
        work.dir = round.dir,
        lastal.command = lastal.command,
        headers = headers,
        min.match.percent = min.match.percent,
        min.match.length = min.match.length,
        min.match.coverage = min.match.coverage,
        max.extension = max.extension,
        max.target.hits = max.target.hits,
        threads = threads,
        quiet = quiet)

      if (length(round.best) == 0) {
        print(paste0(sample, " round ", round, ": no contig passed the target filters."))
        break
      }

      # A later round must not lose a target that an earlier round recovered
      best.contigs = .keepLonger(best.contigs, round.best)

      print(paste0(sample, " round ", round, ": ", length(round.best),
                   " targets assembled, ", length(best.contigs), " total."))

      # Later rounds only extend. A target that gives no contig in one round also
      # gives no contig in the next round from the same bait.
      bait.set = list(
        seqs = stats::setNames(best.contigs,
                               sprintf("bait%06d", seq_along(best.contigs))),
        table = data.frame(bait = sprintf("bait%06d", seq_along(best.contigs)),
                           locus = names(best.contigs),
                           source = "round",
                           stringsAsFactors = FALSE))

    } # end round loop

    unlink(paste0(sample.dir, c("/all_R1.fastq", "/all_R2.fastq")))

    if (is.null(best.contigs) == TRUE || length(best.contigs) == 0) {
      print(paste0(sample, ": nothing was recovered. The previous assembly is saved unchanged."))
      if (is.null(old.contigs) == FALSE) Biostrings::writeXStringSet(old.contigs, out.file)
      next
    }

    final.contigs = .mergeBinnedAssembly(old.contigs = old.contigs,
                                         new.contigs = best.contigs,
                                         locus.names = locus.names,
                                         multi.copy = multi.copy)

    names(final.contigs) = make.unique(names(final.contigs), sep = "_")
    Biostrings::writeXStringSet(final.contigs, out.file)

    .writeBinnedStats(sample = sample,
                      sample.dir = sample.dir,
                      bait.table = stats.table,
                      locus.names = locus.names,
                      new.contigs = best.contigs,
                      old.contigs = old.contigs)

    print(paste0(sample, " finished: ", length(best.contigs),
                 " targets assembled from bins, ", length(final.contigs),
                 " contigs in the output."))

  }, error = function(e) {
    print(paste0(sample.names[i], " failed: ", conditionMessage(e)))
  })
  } # end sample loop

  unlink(db.dir, recursive = TRUE)

} # end function
