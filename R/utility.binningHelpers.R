# Internal helpers for assembleBinnedTargets. These functions are not exported.
# They find the baits, extract the reads of one bin, assemble it, and merge the
# result into the assembly of the main pipeline.
#
# Every step uses an established program. samtools does the read extraction and
# the FASTQ conversion, LAST does the divergent search, and bwa does the
# mapping. There is no custom parsing of alignment records.


# Finds the read pair files of one sample. Reads are selected by name, not by
# position, so every lane is used and an extra file cannot shift the pair.
.pairSampleReads = function(sample.read.dir = NULL) {

  if (dir.exists(sample.read.dir) == FALSE) return(NULL)

  set.reads = list.files(sample.read.dir, full.names = TRUE)
  set.reads = set.reads[grep("fastq|fq", basename(set.reads))]
  if (length(set.reads) == 0) return(NULL)

  read1 = sort(set.reads[grep("_1\\.f|-1\\.f|_R1[_.-]|-R1[_.-]|READ1", basename(set.reads))])
  read2 = sort(set.reads[grep("_2\\.f|-2\\.f|_R2[_.-]|-R2[_.-]|READ2", basename(set.reads))])

  # fastp writes the merged read of an overlapping pair to READ3. It carries the
  # whole insert, so it holds the flank that the mate of an unmerged pair gives.
  read3 = sort(set.reads[grep("_3\\.f|-3\\.f|_R3[_.-]|-R3[_.-]|READ3", basename(set.reads))])

  if (length(read1) == 0 || length(read1) != length(read2)) return(NULL)

  return(list(read1 = read1, read2 = read2, read3 = read3))
}#end .pairSampleReads


# TRUE when a contig source was given but this sample has no file in any of
# them. This is a file test only, so it can run before the lanes are joined and
# before any read is touched. It catches the common case, a sample that is
# absent from the assembly. `.sampleHasNoContigs` catches the rest, a file that
# exists but yields nothing.
.sampleHasNoContigFile = function(sample = NULL,
                                  assembly.directory = NULL,
                                  draft.assembly.directory = NULL) {

  has.file = function(one.dir) {
    if (is.null(one.dir) == TRUE) return(FALSE)
    file.exists(paste0(sub("/+$", "", one.dir), "/", sample, ".fa"))
  }

  has.source = is.null(assembly.directory) == FALSE ||
               is.null(draft.assembly.directory) == FALSE

  return(has.source == TRUE &&
         has.file(assembly.directory) == FALSE &&
         has.file(draft.assembly.directory) == FALSE)
}#end .sampleHasNoContigFile


# TRUE when a contig source was given but this sample holds nothing in it. Such
# a sample cannot be extended or patched, and every target falls into the rescue
# pool, which is days of cap3 on a large read set for a result that reference
# baits alone would give. A run given no contig source at all is a deliberate
# reference-only run and is left alone.
.sampleHasNoContigs = function(own.contigs = NULL,
                               draft.contigs = NULL,
                               assembly.directory = NULL,
                               draft.assembly.directory = NULL) {

  has.source = is.null(assembly.directory) == FALSE ||
               is.null(draft.assembly.directory) == FALSE

  return(has.source == TRUE &&
         length(own.contigs) == 0 &&
         length(draft.contigs) == 0)
}#end .sampleHasNoContigs


# Strips the suffix that make.unique added to a duplicated target name. A target
# name can hold an underscore itself, so a name is only cut when the cut form is
# a real target.
.baseTargetName = function(contig.names = NULL,
                           locus.names = NULL) {

  base.names = contig.names
  needs.cut  = base.names %in% locus.names == FALSE
  cut.names  = sub("_[0-9]+$", "", base.names[needs.cut])
  cut.names[cut.names %in% locus.names == FALSE] = base.names[needs.cut][cut.names %in% locus.names == FALSE]
  base.names[needs.cut] = cut.names

  return(base.names)
}#end .baseTargetName


# Keeps the longest contig for each target. The input contigs can carry the
# make.unique suffix of the main pipeline, so the names are cut back first.
.longestPerTarget = function(contigs = NULL,
                             locus.names = NULL) {

  if (is.null(contigs) == TRUE || length(contigs) == 0) return(NULL)

  base.names = .baseTargetName(names(contigs), locus.names)
  keep       = base.names %in% locus.names
  contigs    = contigs[keep]
  base.names = base.names[keep]
  if (length(contigs) == 0) return(NULL)

  order.index = order(base.names, -Biostrings::width(contigs))
  contigs     = contigs[order.index]
  base.names  = base.names[order.index]

  contigs = contigs[duplicated(base.names) == FALSE]
  names(contigs) = base.names[duplicated(base.names) == FALSE]

  return(contigs)
}#end .longestPerTarget


# Searches the draft assembly of one sample against the targets with LAST.
#
# This is the step that handles divergence. A target locus assembles in the
# draft assembly whatever its divergence, because assembly needs no reference.
# The locus is then lost at the identification step, where blastn cannot match a
# contig that is 35 percent divergent from the probe. LAST trains a scoring
# matrix on this sample and finds those contigs. Numbers in HANDOFF.md.
.identifyDraftContigs = function(draft.file = NULL,
                                 last.db = NULL,
                                 work.dir = NULL,
                                 lastal.command = NULL,
                                 headers = NULL,
                                 min.match.percent = 60,
                                 min.match.length = 50,
                                 min.match.coverage = 30,
                                 threads = 1,
                                 quiet = TRUE) {

  hits.file = paste0(work.dir, "/last-hits.txt")

  # last-train is not used. It learns one substitution rate for the whole input,
  # and a draft assembly holds mostly near-identical contigs. The trained matrix
  # then rejects the few divergent contigs, which are the ones this step must
  # find. The default scores find both. Numbers in HANDOFF.md.
  #
  # BlastTab+ gives the same first 14 columns as the BLAST format of the main
  # pipeline, so the same filters apply. grep removes the comment header, and it
  # returns 1 when there is no hit, which is not an error here.
  .runCommand(paste0(lastal.command, " -P ", threads,
                     " -f BlastTab+ ", shQuote(last.db), " ", shQuote(draft.file),
                     " | grep -v '^#' > ", shQuote(hits.file), " || true"),
              quiet = quiet, task = "LAST search", keep.stdout = TRUE)

  if (file.exists(hits.file) == FALSE || file.size(hits.file) == 0) {
    unlink(hits.file)
    return(NULL)
  }

  hit.data = data.table::fread(hits.file, sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE)
  unlink(hits.file)
  if (ncol(hit.data) < length(headers)) return(NULL)
  data.table::setnames(hit.data, seq_along(headers), headers)

  hit.data = hit.data[hit.data$matches >= min.match.length, ]
  hit.data = hit.data[hit.data$pident >= min.match.percent, ]
  hit.data = hit.data[hit.data$matches >= ((min.match.coverage / 100) * hit.data$tLen), ]
  if (nrow(hit.data) == 0) return(NULL)

  # One draft contig per target, the strongest hit first
  data.table::setorderv(hit.data, c("tName", "bitscore", "pident", "qLen"),
                        order = c(1L, -1L, -1L, -1L))
  best.data = hit.data[duplicated(hit.data$tName) == FALSE, ]

  draft.contigs = Biostrings::readDNAStringSet(draft.file)
  names(draft.contigs) = gsub(" .*", "", names(draft.contigs))

  keep.index = match(best.data$qName, names(draft.contigs))
  best.data  = best.data[is.na(keep.index) == FALSE, ]
  keep.index = keep.index[is.na(keep.index) == FALSE]
  if (length(keep.index) == 0) return(NULL)

  found.contigs = draft.contigs[keep.index]
  names(found.contigs) = best.data$tName

  return(found.contigs)
}#end .identifyDraftContigs


# Chooses one bait for each target and gives it a plain identifier.
#
# The bait is a sequence of the sample when the sample has one that covers
# enough of the target, because reads then map to their own sample at full
# identity and bwa is sensitive. A short fragment recruits only across itself,
# so the reference is used instead. Bait names are numbered identifiers, because
# samtools reads a region as name:start-end and a target name can hold a colon.
.buildBaitTable = function(reference.seqs = NULL,
                           own.contigs = NULL,
                           draft.contigs = NULL,
                           seed.contigs = NULL,
                           target.names = NULL,
                           bait.source = "hybrid",
                           min.bait.coverage = 0.5) {

  # Three places can hold a sequence of this sample for a target. The longest of
  # them recruits the most reads.
  sample.sets = list(contig = own.contigs,
                     draft  = draft.contigs,
                     rescue = seed.contigs)

  bait.seqs = vector("list", length(target.names))
  bait.tags = character(length(target.names))

  for (i in seq_along(target.names)) {
    target = target.names[i]
    reference.seq = reference.seqs[[target]]

    sample.seq = NULL
    sample.tag = ""
    for (tag in names(sample.sets)) {
      one.set = sample.sets[[tag]]
      if (is.null(one.set) == TRUE) next
      if (target %in% names(one.set) == FALSE) next
      candidate = one.set[[target]]
      if (is.null(sample.seq) == TRUE || length(candidate) > length(sample.seq)) {
        sample.seq = candidate
        sample.tag = tag
      }
    }

    use.sample = is.null(sample.seq) == FALSE &&
                 bait.source != "reference" &&
                 length(sample.seq) >= (min.bait.coverage * length(reference.seq))

    if (use.sample == TRUE) {
      bait.seqs[[i]] = sample.seq
      bait.tags[i]   = sample.tag
    } else {
      bait.seqs[[i]] = reference.seq
      bait.tags[i]   = "reference"
    }
  }

  bait.set = Biostrings::DNAStringSet(bait.seqs)
  names(bait.set) = sprintf("bait%06d", seq_along(target.names))

  bait.table = data.frame(bait = names(bait.set),
                          locus = target.names,
                          source = bait.tags,
                          stringsAsFactors = FALSE)

  return(list(seqs = bait.set, table = bait.table))
}#end .buildBaitTable


# Extracts the reads of one bin and assembles them, then deletes the files.
#
# samtools reads the bin straight from the indexed BAM. A pair whose mate did
# not map is stored at the position of the mapped mate, so the mate comes with
# it and supplies the flanking sequence. The output directory of a small megahit
# run has about 27 files. Tens of thousands of bins therefore use the inode quota
# of a shared filesystem before the disk quota.
.assembleOneBin = function(bait = NULL,
                           bam.file = NULL,
                           bin.dir = NULL,
                           samtools.command = NULL,
                           megahit.command = NULL,
                           kmer.values = NULL,
                           subsample = NULL,
                           memory = 4,
                           min.contig.length = 100,
                           quiet = TRUE) {

  read1        = paste0(bin.dir, "/", bait, "_R1.fastq")
  read2        = paste0(bin.dir, "/", bait, "_R2.fastq")
  reads.single = paste0(bin.dir, "/", bait, "_S.fastq")
  assembly.dir = paste0(bin.dir, "/megahit_", bait)
  empty.result = Biostrings::DNAStringSet()

  view.options = ""
  if (is.null(subsample) == FALSE) view.options = paste0(" -s ", subsample)

  extract.status = suppressWarnings(system(
    paste0(samtools.command, " view -b", view.options, " ",
           shQuote(bam.file), " ", bait,
           " | ", samtools.command, " collate -O -u - ",
           shQuote(paste0(bin.dir, "/collate_", bait)),
           " | ", samtools.command, " fastq -N",
           " -1 ", shQuote(read1), " -2 ", shQuote(read2),
           " -0 /dev/null -s ", shQuote(reads.single), " -"),
    ignore.stdout = TRUE, ignore.stderr = TRUE))

  # A merged read maps on its own, so it lands here as a singleton. Discarding it
  # would drop every overlapping pair of the library.
  has.size = function(f) file.exists(f) == TRUE && file.size(f) > 0
  have.pairs  = has.size(read1) && has.size(read2)
  have.single = has.size(reads.single)

  if (extract.status != 0 || (have.pairs == FALSE && have.single == FALSE)) {
    unlink(c(read1, read2, reads.single))
    return(empty.result)
  }

  read.args = ""
  if (have.pairs == TRUE) {
    read.args = paste0(" -1 ", shQuote(read1), " -2 ", shQuote(read2))
  }
  if (have.single == TRUE) {
    read.args = paste0(read.args, " -r ", shQuote(reads.single))
  }

  # megahit assembles a bin that SPAdes drops. SPAdes fits a k-mer coverage model
  # and returns nothing below about 30 read pairs, which is a third to a half of
  # all bins. megahit also runs flat with bin depth. Numbers in HANDOFF.md.
  # megahit refuses an output directory that exists, so a leftover one is removed.
  # -m is in bytes here. The default is a fraction of total RAM, which is wrong
  # when many bins run at once.
  unlink(assembly.dir, recursive = TRUE)
  assembly.status = suppressWarnings(system(
    paste0(megahit.command,
           " --k-list ", paste(kmer.values, collapse = ","),
           " -t 1 -m ", format(memory * 1e9, scientific = FALSE),
           " --min-contig-len ", min.contig.length,
           read.args,
           " -o ", shQuote(assembly.dir)),
    ignore.stdout = TRUE, ignore.stderr = TRUE))

  contig.file = paste0(assembly.dir, "/final.contigs.fa")
  result = empty.result

  if (assembly.status == 0 && file.exists(contig.file) == TRUE &&
      file.size(contig.file) > 0) {
    result = Biostrings::readDNAStringSet(contig.file)
    result = result[Biostrings::width(result) >= min.contig.length]
    if (length(result) > 0) {
      names(result) = paste0(bait, "_", gsub(" .*", "", names(result)))
    }
  }

  unlink(c(assembly.dir, read1, read2, reads.single), recursive = TRUE)

  return(result)
}#end .assembleOneBin


# Flags a sequence when its most frequent 2-mer covers more than half of it.
# Such a sequence extended out of the target into a simple repeat. It is longer
# than the correct contig, so the length rule would otherwise keep it.
.lowComplexity = function(seqs = NULL,
                          max.dimer.fraction = 0.5) {

  if (length(seqs) == 0) return(logical(0))

  dimer.counts = Biostrings::oligonucleotideFrequency(seqs, width = 2)
  totals       = rowSums(dimer.counts)
  totals[totals == 0] = 1

  return((apply(dimer.counts, 1, max) / totals) > max.dimer.fraction)
}#end .lowComplexity


# Merges two named sequence sets and keeps the longer sequence for each name.
# Names that only one set holds are kept unchanged.
.keepLonger = function(current = NULL,
                       candidate = NULL) {

  if (is.null(current) == TRUE || length(current) == 0) return(candidate)
  if (is.null(candidate) == TRUE || length(candidate) == 0) return(current)

  shared = intersect(names(current), names(candidate))
  if (length(shared) > 0) {
    current.width   = Biostrings::width(current)[match(shared, names(current))]
    candidate.width = Biostrings::width(candidate)[match(shared, names(candidate))]
    replace.names   = shared[candidate.width > current.width]
    current = current[names(current) %in% replace.names == FALSE]
  }

  new.names = names(candidate)[names(candidate) %in% names(current) == FALSE]

  return(append(current, candidate[names(candidate) %in% new.names]))
}#end .keepLonger


# Searches the new contigs against the targets with LAST. Keeps one contig per
# target. A contig is only given to the locus of the reads that assembled it, so
# a repeat from one bin cannot go to a different target.
#
# LAST does this search, not blastn. A contig rescued from a divergent locus is
# as divergent from its probe as the reads were, and blastn cannot match it. The
# locus was recovered and then dropped at this gate. Numbers in HANDOFF.md.
.filterBinnedContigs = function(contigs = NULL,
                                bait.table = NULL,
                                bait.widths = NULL,
                                db.prefix = NULL,
                                work.dir = NULL,
                                lastal.command = NULL,
                                headers = NULL,
                                min.match.percent = 60,
                                min.match.length = 50,
                                min.match.coverage = 30,
                                max.extension = 1000,
                                max.target.hits = 5,
                                threads = 1,
                                quiet = TRUE) {

  empty.result = Biostrings::DNAStringSet()
  if (length(contigs) == 0) return(empty.result)

  query.file = paste0(work.dir, "/binned-contigs.fa")
  blast.out  = paste0(work.dir, "/binned-hits.txt")
  Biostrings::writeXStringSet(contigs, query.file)

  .lastSearch(query.file = query.file,
              db.prefix = db.prefix,
              out.file = blast.out,
              lastal.command = lastal.command,
              threads = threads,
              quiet = quiet)

  unlink(query.file)
  if (file.exists(blast.out) == FALSE || file.size(blast.out) == 0) {
    unlink(blast.out)
    return(empty.result)
  }

  match.data = data.table::fread(blast.out, sep = "\t", header = FALSE,
                                 stringsAsFactors = FALSE)
  unlink(blast.out)
  data.table::setnames(match.data, headers)

  filt.data = match.data[match.data$matches >= min.match.length, ]
  filt.data = filt.data[filt.data$pident >= min.match.percent, ]
  filt.data = filt.data[filt.data$matches >= ((min.match.coverage / 100) * filt.data$tLen), ]
  if (nrow(filt.data) == 0) return(empty.result)

  # A contig that passes against many targets at once is a shared repeat.
  # Distinct targets are counted, not rows: one target often gives several HSPs.
  hit.counts = tapply(filt.data$tName, filt.data$qName,
                      function(x) length(unique(x)))
  busy.names = names(hit.counts)[hit.counts > max.target.hits]
  filt.data  = filt.data[filt.data$qName %in% busy.names == FALSE, ]
  if (nrow(filt.data) == 0) return(empty.result)

  # .assembleOneBin prefixes every contig with its bait name. What follows is the
  # assembler's own name for the contig, which is not the same between
  # assemblers, so only the prefix is matched.
  contig.bait = sub("^(bait[0-9]+)_.*$", "\\1", filt.data$qName)
  filt.data$binLocus = bait.table$locus[match(contig.bait, bait.table$bait)]
  filt.data = filt.data[is.na(filt.data$binLocus) == FALSE, ]
  filt.data = filt.data[filt.data$binLocus == filt.data$tName, ]
  if (nrow(filt.data) == 0) return(empty.result)

  # Bitscore first, then length as the tie break. This is the rule the rest of
  # the package uses. Length alone can prefer a contig that extended into a
  # repeat over a better supported one. setorderv takes column names, which
  # prevents the data.table global variable notes that setorder causes in
  # R CMD check.
  data.table::setorderv(filt.data,
                        c("tName", "bitscore", "qLen", "pident", "evalue"),
                        order = c(1L, -1L, -1L, -1L, 1L))
  best.data = filt.data[duplicated(filt.data$tName) == FALSE, ]

  keep.index = match(best.data$qName, names(contigs))
  best.data  = best.data[is.na(keep.index) == FALSE, ]
  keep.index = keep.index[is.na(keep.index) == FALSE]
  if (length(keep.index) == 0) return(empty.result)

  new.contigs = contigs[keep.index]
  names(new.contigs) = best.data$tName

  # Growth cap. A contig can extend by max.extension on each side of its bait.
  # Without the cap, a round can extend into a transposable element. The length
  # rule would then keep that contig.
  if (is.null(bait.widths) == FALSE) {
    allowed = bait.widths[names(new.contigs)] + (2 * max.extension)
    allowed[is.na(allowed)] = Inf
    new.contigs = new.contigs[Biostrings::width(new.contigs) <= allowed]
  }
  if (length(new.contigs) == 0) return(empty.result)

  new.contigs = new.contigs[.lowComplexity(new.contigs) == FALSE]

  return(new.contigs)
}#end .filterBinnedContigs


# Merges the binned contigs into the assembly of the main pipeline. A target with
# one contig keeps the longer of the two sequences. The default leaves a target
# with more than one contig unchanged. One binned contig cannot represent two
# copies. reduceRedundancy also collapses the near-identical haplotype contigs
# before this step, so most of the other targets are paralogs.
.mergeBinnedAssembly = function(old.contigs = NULL,
                                new.contigs = NULL,
                                locus.names = NULL,
                                multi.copy = "keep") {

  if (is.null(old.contigs) == TRUE || length(old.contigs) == 0) return(new.contigs)

  base.names = .baseTargetName(names(old.contigs), locus.names)
  old.widths = Biostrings::width(old.contigs)
  new.widths = Biostrings::width(new.contigs)

  keep.old = rep(TRUE, length(old.contigs))
  keep.new = rep(FALSE, length(new.contigs))

  for (i in seq_along(new.contigs)) {
    rows = which(base.names == names(new.contigs)[i])

    # The target is new to this sample
    if (length(rows) == 0) {
      keep.new[i] = TRUE
      next
    }

    # The target has more than one copy, so the binned contig cannot replace it
    if (length(rows) > 1 && multi.copy == "keep") next

    longest = rows[which.max(old.widths[rows])]
    if (new.widths[i] > old.widths[longest]) {
      keep.new[i]      = TRUE
      keep.old[longest] = FALSE
    }
  }

  return(append(old.contigs[keep.old], new.contigs[keep.new]))
}#end .mergeBinnedAssembly


# Writes the per-sample report. Use this table to judge a run. It gives the reads
# in each bin, the bait that recruited them, and the length before and after.
.writeBinnedStats = function(sample = NULL,
                             sample.dir = NULL,
                             bait.table = NULL,
                             locus.names = NULL,
                             new.contigs = NULL,
                             old.contigs = NULL) {

  if (is.null(new.contigs) == TRUE || length(new.contigs) == 0) return(invisible(NULL))

  old.width = stats::setNames(rep(0, length(locus.names)), locus.names)
  if (is.null(old.contigs) == FALSE && length(old.contigs) > 0) {
    base.names = .baseTargetName(names(old.contigs), locus.names)
    best.old   = tapply(Biostrings::width(old.contigs), base.names, max)
    shared     = intersect(names(best.old), locus.names)
    old.width[shared] = best.old[shared]
  }

  target = names(new.contigs)
  row.index = match(target, bait.table$locus)

  stats.table = data.frame(
    sample         = sample,
    target         = target,
    baitSource     = bait.table$source[row.index],
    readsBinned    = bait.table$reads[row.index],
    binnedLength   = Biostrings::width(new.contigs),
    previousLength = as.numeric(old.width[target]),
    stringsAsFactors = FALSE
  )
  stats.table$finalLength = pmax(stats.table$binnedLength, stats.table$previousLength)
  stats.table$source = ifelse(stats.table$binnedLength > stats.table$previousLength,
                              "binned", "previous")

  write.table(stats.table,
              file = paste0(sample.dir, "/", sample, "_binned-stats.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  return(invisible(stats.table))
}#end .writeBinnedStats


# Adds one row for a finished sample to a CSV that grows across a batch. The row
# is written as each sample finishes, so a batch that stops early keeps the rows
# it earned. A rerun of one sample replaces its row instead of adding a second.
.appendBinnedSummary = function(sample = NULL,
                                locus.stats = NULL,
                                bait.table = NULL,
                                old.contigs = NULL,
                                final.contigs = NULL,
                                rescue.pool = 0,
                                rescue.seeds = 0,
                                divergent.pool = 0,
                                divergent.recovered = 0,
                                bins.per.round = integer(0),
                                targets.per.round = integer(0),
                                minutes = NA_real_,
                                log.directory = "logs") {

  if (is.null(locus.stats) == TRUE || nrow(locus.stats) == 0) return(invisible(NULL))

  had        = locus.stats$previousLength > 0
  gain       = locus.stats$finalLength - locus.stats$previousLength
  extended   = gain > 0 & had
  bait.count = function(tag) sum(bait.table$source == tag)

  previous.bp = if (is.null(old.contigs)) 0 else sum(Biostrings::width(old.contigs))
  final.bp    = if (is.null(final.contigs)) 0 else sum(Biostrings::width(final.contigs))

  summary.row = data.frame(
    sample                = sample,
    targetsBinned         = nrow(locus.stats),
    contigsOut            = if (is.null(final.contigs)) 0 else length(final.contigs),
    targetsWithContig     = sum(had),
    targetsNew            = sum(had == FALSE),
    targetsExtended       = sum(extended),
    percentExtended       = if (sum(had) > 0) round(100 * sum(extended) / sum(had), 1) else NA_real_,
    medianBpAdded         = if (any(extended)) stats::median(gain[extended]) else 0,
    totalBpExtended       = sum(gain[extended]),
    previousBp            = previous.bp,
    finalBp               = final.bp,
    bpGained              = final.bp - previous.bp,
    baitContig            = bait.count("contig"),
    baitDraft             = bait.count("draft"),
    baitReference         = bait.count("reference"),
    baitRescue            = bait.count("rescue"),
    rescueMissingPool     = rescue.pool,
    rescueMissingSeeds    = rescue.seeds,
    divergentPool         = divergent.pool,
    divergentRecovered    = divergent.recovered,
    rounds                = length(targets.per.round),
    round1Bins            = if (length(bins.per.round) > 0) bins.per.round[1] else NA_integer_,
    round1Targets         = if (length(targets.per.round) > 0) targets.per.round[1] else NA_integer_,
    minutes               = round(minutes, 1),
    stringsAsFactors      = FALSE
  )

  dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
  out.csv = paste0(sub("/+$", "", log.directory), "/assembleBinnedTargets_summary.csv")

  if (file.exists(out.csv) == TRUE) {
    existing = utils::read.csv(out.csv, stringsAsFactors = FALSE)
    if ("sample" %in% colnames(existing) == TRUE) {
      existing = existing[existing$sample %in% summary.row$sample == FALSE, , drop = FALSE]
    }
    # An older file can hold a different set of columns. Fill both sides so the
    # rows of a previous version are kept rather than dropped.
    for (m in setdiff(colnames(existing), colnames(summary.row))) summary.row[[m]] = NA
    for (m in setdiff(colnames(summary.row), colnames(existing))) existing[[m]] = NA
    if (nrow(existing) > 0) {
      summary.row = rbind(existing[, colnames(summary.row), drop = FALSE], summary.row)
    }
  }

  utils::write.csv(summary.row, file = out.csv, row.names = FALSE)

  return(invisible(out.csv))
}#end .appendBinnedSummary


# Recovers targets that have no sequence in this sample at all.
#
# cap3 assembles the rescued reads of one target. The reads are unpaired and few,
# so an overlap assembler fits where a de Bruijn one does not. cap3 writes its
# output beside the input, so each target gets its own file name. The read count
# is capped before cap3 sees it, because cap3 compares every read with every
# other one. Numbers in HANDOFF.md.
.assembleOneRescue = function(target = NULL,
                              bam.file = NULL,
                              target.dir = NULL,
                              samtools.command = NULL,
                              cap3.command = NULL,
                              subsample = NULL) {

  reads.file = paste0(target.dir, "/", target, ".fa")
  empty.result = Biostrings::DNAStringSet()
  cap3.files = paste0(reads.file, c("", ".cap.contigs", ".cap.contigs.links",
                                    ".cap.contigs.qual", ".cap.ace", ".cap.info",
                                    ".cap.singlets"))

  # ignore.stdout must stay FALSE. R adds its own redirect for it, which would
  # win over the one that writes reads.file and leave the file empty.
  view.options = ""
  if (is.null(subsample) == FALSE) view.options = paste0(" -s ", subsample)

  extract.status = suppressWarnings(system(
    paste0(samtools.command, " view -b", view.options, " ", shQuote(bam.file), " ", target,
           " | ", samtools.command, " fasta - > ", shQuote(reads.file)),
    ignore.stdout = FALSE, ignore.stderr = TRUE))

  if (extract.status != 0 || file.exists(reads.file) == FALSE ||
      file.size(reads.file) == 0) {
    unlink(cap3.files)
    return(empty.result)
  }

  cap3.status = suppressWarnings(system(paste0(cap3.command, " ", shQuote(reads.file)),
                                        ignore.stdout = TRUE, ignore.stderr = TRUE))

  contig.file = paste0(reads.file, ".cap.contigs")
  result = empty.result

  if (cap3.status == 0 && file.exists(contig.file) == TRUE &&
      file.size(contig.file) > 0) {
    result = Biostrings::readDNAStringSet(contig.file)
    if (length(result) > 0) {
      names(result) = paste0(target, "_", gsub(" .*", "", names(result)))
    }
  }

  unlink(cap3.files)

  return(result)
}#end .assembleOneRescue


# bwa cannot recruit a read that is 35 percent divergent from its bait. On a test
# of 207 read pairs from such a locus, bwa mapped none and LAST aligned 293 of
# the 414 reads. Numbers in HANDOFF.md. LAST therefore does the recruiting here.
#
# The recruited reads are assembled one target at a time. A full run rescues about
# 24,000 targets and recruits millions of alignments, which no single assembly can
# hold. LAST reports one alignment per read and sets no pair flags, so these reads
# are unpaired and shallow: a median of about 9 per target. cap3 suits that. It
# overlaps reads in pairs and applies no coverage model, so it returns a seed from
# as few as two reads where a de Bruijn assembler returns nothing. Numbers in
# HANDOFF.md. The seed only has to bait the paired rounds that follow, and those
# rounds recover the mates and the flanking sequence.
.rescueMissingTargets = function(missing.seqs = NULL,
                                 read.files = NULL,
                                 work.dir = NULL,
                                 lastdb.command = NULL,
                                 lastal.command = NULL,
                                 mafconvert.command = NULL,
                                 samtools.command = NULL,
                                 cap3.command = NULL,
                                 headers = NULL,
                                 min.match.percent = 60,
                                 min.match.length = 50,
                                 min.match.coverage = 30,
                                 max.reads = 6000,
                                 threads = 1,
                                 quiet = TRUE) {

  if (is.null(missing.seqs) == TRUE || length(missing.seqs) == 0) return(NULL)

  # The working files go under the output directory, never in a system temp
  # place. The MAF of a full read set is several GB, so cleanup runs on every
  # exit path, including a failed external program.
  rescue.dir = paste0(work.dir, "/rescue")
  unlink(rescue.dir, recursive = TRUE)
  dir.create(rescue.dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(rescue.dir, recursive = TRUE), add = TRUE)

  missing.file = paste0(rescue.dir, "/missing-targets.fa")
  Biostrings::writeXStringSet(missing.seqs, missing.file)
  .lastBuildDB(reference.file = missing.file,
               db.prefix = paste0(rescue.dir, "/missing_db"),
               lastdb.command = lastdb.command,
               threads = threads,
               quiet = quiet)

  # The reads are renamed on the way in. Only the sequence matters here, and a
  # plain name cannot clash with the region syntax of samtools.
  reads.fasta = paste0(rescue.dir, "/reads.fa")
  .runCommand(paste0("gzip -cdf ", paste(shQuote(read.files), collapse = " "),
                     " | awk 'NR%4==1{printf \">r%d\\n\", ++n} NR%4==2{print}'",
                     " > ", shQuote(reads.fasta)),
              quiet = quiet, task = "read conversion", keep.stdout = TRUE)

  # LAST writes MAF, maf-convert writes SAM, and samtools collects the reads.
  # No alignment record is parsed here. The MAF goes to a file because
  # maf-convert -d does not accept a pipe. The file is deleted straight after.
  rescue.maf = paste0(rescue.dir, "/rescue.maf")
  rescue.bam = paste0(rescue.dir, "/rescue.bam")

  .runCommand(paste0(samtools.command, " faidx ", shQuote(missing.file)),
              quiet = quiet, task = "target indexing")
  .runCommand(paste0(lastal.command, " -P ", threads, " -f MAF ",
                     shQuote(paste0(rescue.dir, "/missing_db")), " ", shQuote(reads.fasta),
                     " > ", shQuote(rescue.maf)),
              quiet = quiet, task = "LAST read rescue", keep.stdout = TRUE)
  unlink(reads.fasta)

  # A coordinate sort and an index let each target's reads be read on their own
  .runCommand(paste0(mafconvert.command, " -d sam ", shQuote(rescue.maf),
                     " | ", samtools.command, " view -bt ", shQuote(paste0(missing.file, ".fai")), " -",
                     " | ", samtools.command, " sort -@ ", threads,
                     " -T ", shQuote(paste0(rescue.dir, "/sorting")),
                     " -o ", shQuote(rescue.bam)),
              quiet = quiet, task = "rescue read collection")
  unlink(rescue.maf)
  .runCommand(paste0(samtools.command, " index -@ ", threads, " ", shQuote(rescue.bam)),
              quiet = quiet, task = "rescue read indexing")

  # idxstats counts the reads of every target in one pass. cap3 joins two reads
  # that overlap, so one read cannot give a seed.
  idx.file = paste0(rescue.dir, "/idxstats.txt")
  .runCommand(paste0(samtools.command, " idxstats ", shQuote(rescue.bam),
                     " > ", shQuote(idx.file)),
              quiet = quiet, task = "rescue read counting", keep.stdout = TRUE)

  idx.data = utils::read.table(idx.file, sep = "\t", header = FALSE,
                               stringsAsFactors = FALSE)
  colnames(idx.data) = c("target", "length", "mapped", "unmapped")
  run.targets = idx.data$target[idx.data$target != "*" & idx.data$mapped >= 2]
  if (length(run.targets) == 0) return(NULL)

  target.dir = paste0(rescue.dir, "/targets")
  dir.create(target.dir, showWarnings = FALSE)

  # A repeat attracts reads from the whole library. In one run a single target
  # took 281,448 of the 1,183,875 recruited reads, and cap3 did not finish it.
  # The cap follows max.pairs in the binned path, and samtools -s takes the seed
  # in the integer part and the fraction to keep after the point.
  target.reads = stats::setNames(idx.data$mapped, idx.data$target)

  contig.list = parallel::mclapply(run.targets, function(target) {
    subsample = NULL
    if (max.reads > 0 && target.reads[[target]] > max.reads) {
      subsample = format(42 + (max.reads / target.reads[[target]]), nsmall = 4)
    }
    .assembleOneRescue(target = target,
                       bam.file = rescue.bam,
                       target.dir = target.dir,
                       samtools.command = samtools.command,
                       cap3.command = cap3.command,
                       subsample = subsample)
  }, mc.cores = min(threads, length(run.targets)), mc.preschedule = FALSE)

  contig.list = contig.list[vapply(contig.list, function(x)
    inherits(x, "DNAStringSet") && length(x) > 0, logical(1))]
  if (length(contig.list) == 0) return(NULL)

  contig.file = paste0(rescue.dir, "/rescue-contigs.fa")
  Biostrings::writeXStringSet(do.call(c, unname(contig.list)), contig.file)

  # A seed is filtered on min.match.length, not min.contig.length. cap3 seeds are
  # short, and a seed shorter than the smallest match worth accepting cannot bait
  # anything. min.contig.length stays with the binned contigs, where the coverage
  # rule in .filterBinnedContigs already sets the real floor. Numbers in HANDOFF.md.
  # The same LAST search that assigns a draft contig assigns a rescued contig
  seed.contigs = .identifyDraftContigs(draft.file = contig.file,
                                       last.db = paste0(rescue.dir, "/missing_db"),
                                       work.dir = rescue.dir,
                                       lastal.command = lastal.command,
                                       headers = headers,
                                       min.match.percent = min.match.percent,
                                       min.match.length = min.match.length,
                                       min.match.coverage = min.match.coverage,
                                       threads = threads,
                                       quiet = quiet)

  if (is.null(seed.contigs) == FALSE) {
    seed.contigs = seed.contigs[Biostrings::width(seed.contigs) >= min.match.length]
    if (length(seed.contigs) == 0) seed.contigs = NULL
  }

  return(seed.contigs)
}#end .rescueMissingTargets
