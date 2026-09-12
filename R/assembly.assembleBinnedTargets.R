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
#'   A sample needs an assembly of its own. One that is absent from both
#'   \code{assembly.directory} and \code{draft.assembly.directory} is skipped
#'   with a warning, because every target would fall into the rescue pool and the
#'   step would run for days to give what reference baits alone give. Assemble
#'   such a sample with \code{assembleSpades} first, then bin it.
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
#'   holds the per-sample read folders. The reads must be paired. A directory
#'   that also holds merged reads, which \code{mergePairedEndReads} writes as
#'   READ3, is used in full: the pairs and the merged reads are mapped in
#'   separate passes into one BAM.
#'
#'   Prefer unmerged reads here. Merging costs targets and returns nothing: it
#'   recovered 213 fewer targets on a test sample and gave the same contig on
#'   the targets both found. A capture insert straddles the target edge, so
#'   merging turns an exonic mate and an intronic mate into one half-off-target
#'   query, and a local alignment can only score over the on-target part. The
#'   pair also gets bwa's mate rescue, which a single-end merged read does not.
#'   Merged reads run about 30 percent faster. Numbers in HANDOFF.md.
#'   Default: \code{"decontaminated-reads"}.
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
#' @param log.directory path to the directory for the cross-sample summary CSV.
#'   Default: \code{"logs"}.
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
#'   assembly also lost cannot be recovered.
#'
#'   Setting it makes those targets take the strict path. The seed becomes a
#'   bait, and the target must then gate a bin at \code{min.pairs}, assemble
#'   under megahit, and clear \code{min.contig.length}. It also removes them
#'   from \code{rescue.failed.divergent}, which writes its contig straight to
#'   the output under the match filters alone. On the test sample the strict
#'   path returned 379 targets and the permissive path 1,879, the extra ones a
#'   median of 135 bp. \code{FALSE} therefore recovers many more targets, and
#'   they are short. Numbers in HANDOFF.md. Default: \code{FALSE}.
#'
#' @param rescue.failed.divergent logical. \code{TRUE} recruits reads with
#'   LAST for the targets that no bin produced, after round 1. A bin fails when
#'   bwa recruited fewer than \code{min.pairs} pairs, or when it assembled and
#'   no contig passed the filters. bwa needs about 90 percent identity to
#'   recruit, so a divergent target recruits nothing however deep it is. The
#'   targets that \code{rescue.missing} already searched are not searched again.
#'   The step costs one more pass over the reads, about 8 minutes for a
#'   5 million read sample, whatever the number of targets. It recovered 916 and
#'   656 targets on the two test runs, about 9 percent of the output, which is
#'   the largest single gain in this function. Numbers in HANDOFF.md.
#'   Default: \code{TRUE}.
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
#'   assembly time and gives the assembler no more information. It also caps the
#'   reads of one rescue target, at twice this value, because those reads are
#'   unpaired. A repeat can otherwise take most of the recruited reads and stall
#'   cap3. \code{0} removes both limits. Default: \code{3000}.
#'
#' @param min.contig.length minimum length in base pairs of a binned contig. The
#'   rescue seeds use \code{min.match.length} instead, because cap3 seeds are
#'   short. Default: \code{100}.
#'
#' @param min.match.percent minimum BLAST percent identity of a new contig
#'   against its target. Default: \code{60}.
#'
#' @param min.match.length minimum alignment length in base pairs. It is also
#'   the length floor of the rescue seeds. Default: \code{50}.
#'
#' @param min.match.coverage minimum percentage of the target length that the
#'   match must cover. Default: \code{30}.
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
#'   assembly step runs this many single-threaded megahit jobs at once. When
#'   \code{parallel.samples > 1} each sample receives
#'   \code{floor(threads / parallel.samples)}. Default: \code{1}.
#'
#' @param parallel.samples number of samples to assemble at the same time.
#'   \code{threads} and \code{memory} are divided between them, the way
#'   \code{assembleSpades} divides them. Unlike SPAdes, one sample here already
#'   uses every thread it is given: the bin assembly is one single-threaded
#'   megahit per bin over \code{mc.cores = threads}, some 10,000 of them, and
#'   the cap3 rescue is the same shape. On a test sample that phase was 88 of
#'   101 minutes, so about 87 percent of a run scales with cores on its own.
#'
#'   Raise this to fill a node across a batch, not to make one sample faster.
#'   With 87 percent parallel, one sample on 48 threads runs about 3.6 times
#'   faster than on 8 but costs 22.4 core-hours against 13.5, because the serial
#'   fraction, bwa mem and the samtools sort and the LAST database, holds the
#'   rest idle. Several samples at 8 to 16 threads each keep those cores busy
#'   instead. The per-bin memory does not change when this is raised, because
#'   \code{memory} and \code{threads} are divided together. Default: \code{1}.
#'   Numbers in HANDOFF.md.
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
#' @param cap3.path path to the directory that holds \code{cap3}. Only
#'   \code{rescue.missing = TRUE} and \code{rescue.failed.divergent = TRUE}
#'   need it. Default: \code{NULL}.
#'
#' @param last.path path to the directory that holds \code{lastdb},
#'   \code{lastal} and \code{maf-convert}. Default: \code{NULL}.
#'
#' @param overwrite logical. \code{TRUE} runs every sample again. \code{FALSE}
#'   skips a sample only when its final binned FASTA in \code{binned.directory}
#'   exists and is not empty; any other sample is assembled again from a clean
#'   working directory. Default: \code{FALSE}.
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
#'   It also adds one row per sample to
#'   \code{log.directory/assembleBinnedTargets_summary.csv}. The row is written
#'   as each sample finishes, so a batch that stops early keeps the rows it
#'   earned, and a rerun of one sample replaces its row. The columns cover the
#'   targets recovered, the targets extended and the base pairs added, the bait
#'   sources, both rescue steps, and the run time. Use it to compare samples
#'   across a large batch.
#'
#'   The row is rewritten by reading the whole file and writing it back, so
#'   concurrent jobs must not share a \code{log.directory}. Give each array task
#'   its own and join the files afterwards.
#'
#'   \code{percentExtended} is measured against whatever sat in
#'   \code{assembly.directory}. A trimmed input set raises it without the step
#'   behaving differently, so it is comparable only between runs that began from
#'   the same contigs. \code{medianPreviousLength} and
#'   \code{medianBinnedLength} sit beside it so the baseline is visible.
#'   \code{baitContig} and the other bait counts are what was offered;
#'   \code{targetsFromContig} and its three partners are what came back.
#'   \code{targetsUnderMinLength} counts the contigs that reached the assembly
#'   below \code{min.contig.length}, which the divergent rescue never tests.
#'
#' @export

assembleBinnedTargets = function(read.directory = NULL,
                                 mapping.reads = "decontaminated-reads",
                                 target.markers = NULL,
                                 assembly.directory = NULL,
                                 draft.assembly.directory = NULL,
                                 output.directory = "binned-target-assembly",
                                 binned.directory = NULL,
                                 log.directory = "logs",
                                 locus.set = c("all", "missing"),
                                 bait.source = c("hybrid", "reference"),
                                 min.bait.coverage = 0.5,
                                 rescue.missing = FALSE,
                                 rescue.failed.divergent = TRUE,
                                 iterations = 1,
                                 min.pairs = 6,
                                 max.pairs = 3000,
                                 min.contig.length = 100,
                                 min.match.percent = 60,
                                 min.match.length = 50,
                                 min.match.coverage = 30,
                                 max.extension = 1000,
                                 max.target.hits = 5,
                                 multi.copy = c("keep", "longest"),
                                 kmer.values = c(21, 33, 55, 77, 99),
                                 memory = 8,
                                 threads = 1,
                                 parallel.samples = 1,
                                 bwa.path = NULL,
                                 samtools.path = NULL,
                                 megahit.path = NULL,
                                 cap3.path = NULL,
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

  # Only the two LAST rescue steps use cap3, so they are what needs it
  use.rescue = rescue.missing == TRUE && bait.source != "reference"
  cap3.command = NULL
  if (use.rescue == TRUE || rescue.failed.divergent == TRUE) {
    cap3.command = .toolCommand("cap3", cap3.path)
  }

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

  # Divides the resources between the samples that run at the same time, the way
  # assembleSpades does. The per-bin allocation is unchanged by the split:
  # each megahit still gets floor(memory / threads), because both are divided by
  # parallel.samples. What changes is how many samples are in flight.
  if (!is.numeric(threads) || length(threads) != 1 || !is.finite(threads) || threads < 1 ||
      !is.numeric(memory) || length(memory) != 1 || !is.finite(memory) || memory <= 0 ||
      !is.numeric(parallel.samples) || length(parallel.samples) != 1 ||
      !is.finite(parallel.samples) || parallel.samples < 1) {
    stop("threads, memory, and parallel.samples must be positive finite values.")
  }
  threads = floor(threads)
  parallel.samples = min(floor(parallel.samples), length(sample.names), threads,
                         max(1, floor(memory / 2)))
  thread.cl = floor(threads / parallel.samples)
  mem.cl    = memory / parallel.samples

  if (parallel.samples > 1) {
    print(paste0("Assembling ", parallel.samples, " samples at a time with ",
                 thread.cl, " threads and ", mem.cl, "GB each."))
  }

  # threads and memory are arguments, so they shadow the totals for everything
  # below and every inner step takes this sample's share without further change.
  assemble.one = function(i, threads, memory) {

  assembly.memory = memory / max(1, threads)
  tryCatch({

    sample       = sample.names[i]
    sample.start = Sys.time()
    sample.dir   = paste0(output.directory, "/", sample)
    out.file   = paste0(binned.directory, "/", sample, ".fa")

    # overwrite = FALSE skips a sample only when its binned FASTA exists and is
    # not empty. The working directory is created below, once the sample is known
    # to be assembled.
    if (overwrite == FALSE && file.exists(out.file) == TRUE &&
        file.info(out.file)$size > 0) {
      print(paste0(sample, " already finished, skipping."))
      return(invisible(NULL))
    }

    read.pair = .pairSampleReads(paste0(actual.read.dir, "/", sample))
    if (is.null(read.pair) == TRUE) {
      warning(sample, ": paired reads were not found in ", mapping.reads, ". Skipping.")
      return(invisible(NULL))
    }

    # A sample absent from every contig source cannot be extended or patched, and
    # every target would fall into the rescue pool, which is days of cap3 on a
    # large read set. This is a file test, so the skip costs nothing and happens
    # before the lanes are joined. Numbers in HANDOFF.md.
    if (.sampleHasNoContigFile(sample = sample,
                               assembly.directory = assembly.directory,
                               draft.assembly.directory = draft.assembly.directory) == TRUE) {
      warning(sample, ": no contigs of its own in the assembly or the draft. ",
              "Skipping. Assemble this sample before binning it.")
      return(invisible(NULL))
    }

    # Give the sample a clean working directory before assembly.
    unlink(sample.dir, recursive = TRUE)
    dir.create(sample.dir, recursive = TRUE, showWarnings = FALSE)

    # bwa mem takes one file per mate, so several lanes are joined first
    read1 = read.pair$read1
    read2 = read.pair$read2
    read3 = read.pair$read3
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
    if (length(read3) > 1) {
      joined.read3 = paste0(sample.dir, "/all_R3.fastq")
      .runCommand(paste0("gzip -cdf ", paste(shQuote(read3), collapse = " "),
                         " > ", shQuote(joined.read3)),
                  quiet = quiet, task = "read joining", keep.stdout = TRUE)
      read3 = joined.read3
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

    # A sample missing from the assembly puts every target in the rescue pool.
    # That is days of work for a result reference baits alone would give, so it
    # is skipped before the rescue runs. Numbers in HANDOFF.md.
    if (.sampleHasNoContigs(own.contigs = own.contigs,
                            draft.contigs = draft.contigs,
                            assembly.directory = assembly.directory,
                            draft.assembly.directory = draft.assembly.directory) == TRUE) {
      warning(sample, ": no contigs of its own in the assembly or the draft. ",
              "Skipping. Assemble this sample before binning it.")
      return(invisible(NULL))
    }

    target.names = locus.names
    if (locus.set == "missing" && is.null(own.contigs) == FALSE) {
      target.names = locus.names[locus.names %in% names(own.contigs) == FALSE]
    }
    if (length(target.names) == 0) {
      print(paste0(sample, ": no target needs binning. The previous assembly is saved unchanged."))
      if (is.null(old.contigs) == FALSE) .writeAtomicFasta(old.contigs, out.file)
      return(invisible(NULL))
    }

    # Targets with no sequence in this sample at all. bwa cannot recruit for
    # these when the sample is divergent, so LAST recruits the reads instead.
    seed.contigs = NULL
    rescue.names = character(0)
    if (use.rescue == TRUE) {
      have.names = c(names(own.contigs), names(draft.contigs))
      rescue.names = target.names[target.names %in% have.names == FALSE]

      if (length(rescue.names) > 0) {
        seed.contigs = .rescueMissingTargets(
          missing.seqs = reference.seqs[rescue.names],
          read.files = c(read1, read2, read3),
          work.dir = sample.dir,
          lastdb.command = lastdb.command,
          lastal.command = lastal.command,
          mafconvert.command = mafconvert.command,
          samtools.command = samtools.command,
          cap3.command = cap3.command,
          headers = headers,
          min.match.percent = min.match.percent,
          min.match.length = min.match.length,
          min.match.coverage = min.match.coverage,
          max.reads = max.pairs * 2,
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

    # Counters for the one row this sample adds to the summary CSV
    bins.per.round      = integer(0)
    targets.per.round   = integer(0)
    divergent.pool      = 0
    divergent.recovered = 0

    run.divergent.rescue = function(round) {
      if (rescue.failed.divergent == FALSE || round != 1) return(FALSE)
      recovered.names = if (is.null(best.contigs)) character(0) else names(best.contigs)
      failed.names = target.names[target.names %in% recovered.names == FALSE &
                                  target.names %in% rescue.names == FALSE]
      if (length(failed.names) == 0) return(FALSE)

      failed.contigs = .rescueMissingTargets(
        missing.seqs = reference.seqs[failed.names],
        read.files = c(read1, read2, read3),
        work.dir = sample.dir,
        lastdb.command = lastdb.command,
        lastal.command = lastal.command,
        mafconvert.command = mafconvert.command,
        samtools.command = samtools.command,
        cap3.command = cap3.command,
        headers = headers,
        min.match.percent = min.match.percent,
        min.match.length = min.match.length,
        min.match.coverage = min.match.coverage,
        max.reads = max.pairs * 2,
        threads = threads,
        quiet = quiet)

      best.contigs <<- .keepLonger(best.contigs, failed.contigs)
      divergent.pool <<- length(failed.names)
      divergent.recovered <<- if (is.null(failed.contigs)) 0 else length(failed.contigs)
      print(paste0(sample, ": LAST recovered ", divergent.recovered,
                   " of ", length(failed.names),
                   " targets that no bin produced."))
      divergent.recovered > 0
    }

    update.round.baits = function() {
      bait.set <<- list(
        seqs = stats::setNames(best.contigs,
                               sprintf("bait%06d", seq_along(best.contigs))),
        table = data.frame(bait = sprintf("bait%06d", seq_along(best.contigs)),
                           locus = names(best.contigs), source = "round",
                           stringsAsFactors = FALSE))
    }

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
      # One bwa call takes either two mate files or one single-end file, so the
      # merged reads need a second pass. samtools cat joins the two unsorted BAMs,
      # which share a header because both used this bait index.
      map.bams      = paste0(round.dir, "/map-pe.bam")
      merged.counts = NULL
      .runCommand(paste0(bwa.command, " mem -t ", threads,
                         " -k 15 -B 3 -O 5 -T 25 ",
                         shQuote(bait.file), " ", shQuote(read1), " ", shQuote(read2),
                         " | ", samtools.command, " view -b -o ", shQuote(map.bams), " -"),
                  quiet = quiet, task = "read mapping")

      if (length(read3) > 0) {
        merged.bam = paste0(round.dir, "/map-se.bam")
        .runCommand(paste0(bwa.command, " mem -t ", threads,
                           " -k 15 -B 3 -O 5 -T 25 ",
                           shQuote(bait.file), " ", shQuote(read3),
                           " | ", samtools.command, " view -b -o ", shQuote(merged.bam), " -"),
                    quiet = quiet, task = "merged read mapping")
        map.bams = c(map.bams, merged.bam)

        # A merged read is one record that spans the whole insert, so the pair
        # count below must score it as a pair and not as half of one. The count
        # is taken here, while the merged reads are still in a file of their own.
        merged.file = paste0(round.dir, "/merged-counts.txt")
        .runCommand(paste0(samtools.command, " view -F 4 ", shQuote(merged.bam),
                           " | awk '{ c[$3]++ } END { for (k in c) print k \"\\t\" c[k] }'",
                           " > ", shQuote(merged.file)),
                    quiet = quiet, task = "merged read counting", keep.stdout = TRUE)

        if (file.exists(merged.file) == TRUE && file.info(merged.file)$size > 0) {
          merged.data   = utils::read.table(merged.file, sep = "\t", header = FALSE,
                                            stringsAsFactors = FALSE)
          merged.counts = stats::setNames(as.numeric(merged.data$V2), merged.data$V1)
        }
      }

      .runCommand(paste0(samtools.command, " cat ",
                         paste(shQuote(map.bams), collapse = " "),
                         " | ", samtools.command, " sort -@ ", threads,
                         " -T ", shQuote(paste0(round.dir, "/sorting")),
                         " -o ", shQuote(bam.file)),
                  quiet = quiet, task = "read sorting")
      unlink(map.bams)
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

      # idxstats counts records, and the test below divides by two to get pairs.
      # A merged read is a whole insert in one record, so its record is added a
      # second time. Without this a merged library loses about an eighth of its
      # inserts at the gate. Numbers in HANDOFF.md.
      if (is.null(merged.counts) == FALSE) {
        extra.reads = merged.counts[names(bait.reads)]
        extra.reads[is.na(extra.reads)] = 0
        bait.reads = bait.reads + extra.reads
      }
      round.reads = as.numeric(bait.reads[bait.set$table$bait])
      round.reads[is.na(round.reads)] = 0
      if (round == 1) stats.table$reads = round.reads[match(stats.table$locus, bait.set$table$locus)]

      run.index = which(round.reads >= (min.pairs * 2))
      bins.per.round = c(bins.per.round, length(run.index))
      print(paste0(sample, " round ", round, ": ", length(run.index), " of ",
                   nrow(bait.set$table), " bins hold at least ", min.pairs,
                   " read pairs. Assembling."))

      if (length(run.index) == 0) {
        rescued = run.divergent.rescue(round)
        if (rescued && round < iterations) {
          update.round.baits()
          next
        }
        break
      }

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
        rescued = run.divergent.rescue(round)
        if (rescued && round < iterations) {
          update.round.baits()
          next
        }
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
        rescued = run.divergent.rescue(round)
        if (rescued && round < iterations) {
          update.round.baits()
          next
        }
        break
      }

      # A later round must not lose a target that an earlier round recovered
      best.contigs = .keepLonger(best.contigs, round.best)
      targets.per.round = c(targets.per.round, length(round.best))

      print(paste0(sample, " round ", round, ": ", length(round.best),
                   " targets assembled, ", length(best.contigs), " total."))

      # Targets that no bin produced. bwa needs about 90 percent identity to
      # recruit, so a divergent target recruits nothing however deep it is, and
      # LAST recruits those reads instead. This runs after round 1, so a later
      # round can extend what it finds. The targets rescue.missing already
      # searched are skipped, because the same reads give the same answer.
      run.divergent.rescue(round)

      # Later rounds only extend. A target that gives no contig in one round also
      # gives no contig in the next round from the same bait.
      update.round.baits()

    } # end round loop

    unlink(paste0(sample.dir, c("/all_R1.fastq", "/all_R2.fastq", "/all_R3.fastq")))

    if (is.null(best.contigs) == TRUE || length(best.contigs) == 0) {
      print(paste0(sample, ": nothing was recovered. The previous assembly is saved unchanged."))
      if (is.null(old.contigs) == FALSE) .writeAtomicFasta(old.contigs, out.file)
      return(invisible(NULL))
    }

    final.contigs = .mergeBinnedAssembly(old.contigs = old.contigs,
                                         new.contigs = best.contigs,
                                         locus.names = locus.names,
                                         multi.copy = multi.copy)

    names(final.contigs) = make.unique(names(final.contigs), sep = "_")
    .writeAtomicFasta(final.contigs, out.file)

    locus.stats = .writeBinnedStats(sample = sample,
                                    sample.dir = sample.dir,
                                    bait.table = stats.table,
                                    locus.names = locus.names,
                                    new.contigs = best.contigs,
                                    old.contigs = old.contigs)

    .appendBinnedSummary(sample = sample,
                         locus.stats = locus.stats,
                         bait.table = stats.table,
                         old.contigs = old.contigs,
                         final.contigs = final.contigs,
                         rescue.pool = length(rescue.names),
                         rescue.seeds = if (is.null(seed.contigs)) 0 else length(seed.contigs),
                         divergent.pool = divergent.pool,
                         divergent.recovered = divergent.recovered,
                         bins.per.round = bins.per.round,
                         targets.per.round = targets.per.round,
                         minutes = as.numeric(difftime(Sys.time(), sample.start, units = "mins")),
                         min.contig.length = min.contig.length,
                         log.directory = log.directory)

    print(paste0(sample, " finished: ", length(best.contigs),
                 " targets assembled from bins, ", length(final.contigs),
                 " contigs in the output."))

  }, error = function(e) {
    print(paste0(sample.names[i], " failed: ", conditionMessage(e)))
  })
  }# end assemble.one

  # A warning raised in a forked child never reaches the parent, so each sample
  # reports itself as it finishes rather than through the return value.
  parallel::mclapply(seq_along(sample.names),
                     function(i) assemble.one(i, thread.cl, mem.cl),
                     mc.cores = parallel.samples)

  # Each sample wrote its own row. Join them once, in the parent, where nothing
  # else is writing.
  .mergeBinnedSummaries(log.directory)

  unlink(db.dir, recursive = TRUE)

} # end function
