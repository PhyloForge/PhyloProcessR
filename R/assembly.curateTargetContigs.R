#' @title curateTargetContigs
#'
#' @description Curates the contigs of each sample against the target markers.
#'   The function reduces redundancy with cd-hit-est, matches the targets to the
#'   contigs with LAST, joins the fragments of one target that sit on separate
#'   contigs, and cuts apart a contig that spans more than one target. It writes
#'   one sequence per target for each sample, named after the target.
#'
#' @details This is the structural half of \code{annotateTargets}. It runs at the
#'   end of workflow 2, before variant calling, because the variant caller maps
#'   the reads back to these contigs. A contig that spans two targets collects
#'   the reads of both loci and gives wrong genotypes, and no later step can
#'   repair that. \code{annotateTargets} keeps the other half of the work: the
#'   paralog policy, the sample naming, and the alignment output.
#'
#'   Run this step with \code{retain.paralogs} in mind. This function does not
#'   drop a paralog copy, because the reads of a dropped copy map to the copy
#'   that remains and make false heterozygosity.
#'
#' @param assembly.directory path to a directory of per-sample contig files, one
#'   \code{<sample>.fa} file per sample.
#'
#' @param target.file path to the FASTA file of target markers.
#'
#' @param output.directory path to the directory for the curated contigs.
#'   Default: \code{"curated-contigs"}.
#'
#' @param min.match.percent minimum percent identity of a match. Default:
#'   \code{60}.
#'
#' @param min.match.length minimum alignment length in base pairs. The N padding
#'   that joins two fragments of one target is not counted. Default: \code{50}.
#'
#' @param min.match.coverage minimum percentage of the target length that the
#'   matches must cover. The hits of one target are summed, and the test runs
#'   after the fragments are joined. A target split across two contigs therefore
#'   passes when the two fragments together cover enough of it. The default is
#'   permissive, because a later step can still remove a short locus, but this
#'   step cannot recover one it dropped. The rescue seeds of
#'   \code{assembleBinnedTargets} are short by design, and later rounds build
#'   them out. Default: \code{30}.
#'
#' @param similarity sequence-identity threshold used by \code{cd-hit-est} to
#'   remove redundant contigs before target matching. Lower values collapse
#'   more similar contigs and can merge recent paralogous copies. Must be from
#'   \code{0.8} through \code{1}. Default: \code{0.9}.
#'
#' @param search.method which program matches the target markers to the contigs.
#'   \code{"last"} (default) uses LAST, which matches a contig that is up to
#'   about 35 percent divergent from its target. \code{"blast"} uses
#'   \code{blastn dc-megablast}, which loses a contig past about 25 percent
#'   divergence. Numbers in HANDOFF.md.
#'
#' @param threads number of samples to process at the same time. Default:
#'   \code{1}.
#'
#' @param memory RAM in GB shared by the parallel samples. Default: \code{1}.
#'
#' @param blast.path path to the directory containing the BLAST executables.
#'   Only needed when \code{search.method = "blast"}. Default: \code{NULL}.
#'
#' @param last.path path to the directory containing \code{lastdb} and
#'   \code{lastal}. Only needed when \code{search.method = "last"}. Default:
#'   \code{NULL}.
#'
#' @param cdhit.path path to the directory containing \code{cd-hit-est}.
#'   Default: \code{NULL}.
#'
#' @param overwrite logical. \code{TRUE} runs a sample again when its output
#'   exists. Default: \code{FALSE}.
#'
#' @param quiet logical. \code{TRUE} hides the output of the external programs.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Writes one \code{<sample>.fa} file per
#'   sample to \code{output.directory}, with each sequence named after its
#'   target, and a match table per sample to
#'   \code{logs/sample_logs/<sample>_curation-matches.csv}.
#'
#' @export

curateTargetContigs = function(assembly.directory = NULL,
                               target.file = NULL,
                               output.directory = "curated-contigs",
                               min.match.percent = 60,
                               min.match.length = 50,
                               min.match.coverage = 30,
                               similarity = 0.9,
                               search.method = c("last", "blast"),
                               threads = 1,
                               memory = 1,
                               blast.path = NULL,
                               last.path = NULL,
                               cdhit.path = NULL,
                               overwrite = FALSE,
                               quiet = TRUE) {

  search.method = match.arg(search.method)

  #Initial checks
  if (is.null(assembly.directory) == TRUE) { stop("Please provide a contig directory.") }
  if (dir.exists(assembly.directory) == FALSE) { stop("Contig directory not found.") }
  if (is.null(target.file) == TRUE) { stop("Please provide a target file.") }
  if (file.exists(target.file) == FALSE) { stop("Target file not found.") }
  if (!is.numeric(similarity) || length(similarity) != 1 ||
      !is.finite(similarity) || similarity < 0.8 || similarity > 1) {
    stop("similarity must be a finite value from 0.8 through 1.")
  }

  if (similarity >= 0.95) {
    word.length = 10
  } else if (similarity >= 0.90) {
    word.length = 8
  } else if (similarity >= 0.88) {
    word.length = 7
  } else if (similarity >= 0.85) {
    word.length = 6
  } else {
    word.length = 5
  }

  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists("logs/sample_logs")) {
    dir.create("logs/sample_logs", recursive = TRUE, showWarnings = FALSE)
  }

  fasta.pattern = "\\.(fa|fas|fasta|fna)$"
  file.names = list.files(assembly.directory, pattern = fasta.pattern,
                          ignore.case = TRUE)
  if (length(file.names) == 0) stop("No FASTA files found in the contig directory.")

  if (!is.numeric(threads) || length(threads) != 1 || !is.finite(threads) || threads < 1 ||
      !is.numeric(memory) || length(memory) != 1 || !is.finite(memory) || memory <= 0) {
    stop("threads and memory must be positive finite values.")
  }
  threads = min(floor(threads), length(file.names))
  cdhit.command = .toolCommand("cd-hit-est", cdhit.path)
  if (search.method == "last") {
    lastdb.command = .toolCommand("lastdb", last.path)
    lastal.command = .toolCommand("lastal", last.path)
  } else {
    makeblastdb.command = .toolCommand("makeblastdb", blast.path)
    blastn.command = .toolCommand("blastn", blast.path)
  }

  #headers for the search results
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  mem.cl = max(1, floor((memory * 1000) / threads))

  results = parallel::mclapply(seq_along(file.names), function(i) {
  tryCatch({

    #Sets up working directories for each species
    sample = sub(fasta.pattern, "", file.names[i], ignore.case = TRUE)

    #Checks if this has been done already (before creating any directory)
    if (overwrite == FALSE){
      existing.file = paste0(output.directory, "/", sample, ".fa")
      if (file.exists(existing.file) == TRUE && file.size(existing.file) > 0) {
        print(paste0(sample, " already finished, skipping. Set overwrite = TRUE to redo."))
        return(NULL)
      }
    }#end
    if (overwrite == TRUE) unlink(paste0(output.directory, "/", sample, ".fa"))

    # Temporary working directory kept inside logs so output.directory stays clean
    species.dir = paste0("logs/sample_logs/", sample)
    if (!dir.exists(species.dir)){ dir.create(species.dir, recursive = TRUE, showWarnings = FALSE) }

    #########################################################################
    # Part A: reduce redundancy
    #########################################################################

    red.file = paste0(species.dir, "/", sample, "_red.fa")
    rename.file = paste0(species.dir, "/", sample, "_rename.fa")
    search.file = paste0(species.dir, "/", sample, "_target-blast-match.txt")
    .runCommand(paste0(cdhit.command, " -i ",
                       shQuote(file.path(assembly.directory, file.names[i])),
                       " -o ", shQuote(red.file), " -p 0 -T 1",
                       " -n ", word.length, " -c ", similarity, " -M ", mem.cl),
                quiet = quiet, task = paste(sample, "cd-hit-est"))
    if (!file.exists(red.file) || file.size(red.file) == 0) {
      stop("cd-hit-est produced no sequences for ", sample, ".")
    }

    ### Read in data
    all.data = Biostrings::readDNAStringSet(file = red.file, format = "fasta")

    names(all.data) = paste0("contig_", seq(seq_along(all.data)))

    # Writes the final loci
    final.loci = as.list(as.character(all.data))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      rename.file,
      nbchar = 1000000, as.string = TRUE, open = "w"
    )

    #########################################################################
    #Part B: Blasting
    #########################################################################

    # The database holds the contigs of this sample and the query is the target
    # file, so a hit names the target first and the contig second.
    # One thread per search, because the samples already run in parallel.
    if (search.method == "last") {
      .lastBuildDB(reference.file = rename.file,
                   db.prefix = paste0(species.dir, "/", sample, "_last_db"),
                   lastdb.command = lastdb.command,
                   threads = 1,
                   quiet = quiet)

      .lastSearch(query.file = target.file,
                  db.prefix = paste0(species.dir, "/", sample, "_last_db"),
                  out.file = search.file,
                  lastal.command = lastal.command,
                  threads = 1,
                  quiet = quiet)
    } else {
      .runCommand(paste0(makeblastdb.command, " -in ", shQuote(rename.file),
        " -parse_seqids -dbtype nucl -out ",
        shQuote(paste0(species.dir, "/", sample, "_nucl-blast_db"))),
        quiet = quiet, task = paste(sample, "BLAST database"))

      .runCommand(paste0(blastn.command, " -task dc-megablast -db ",
        shQuote(paste0(species.dir, "/", sample, "_nucl-blast_db")), " -evalue 0.001",
        " -query ", shQuote(target.file), " -out ", shQuote(search.file),
        " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
        " -num_threads 1"
      ), quiet = quiet, task = paste(sample, "BLAST search"))
    }

    # Remove the search database and the large intermediate contig files
    cleanup.files = list.files(species.dir, full.names = TRUE)
    cleanup.base = basename(cleanup.files)
    cleanup.files = cleanup.files[
      startsWith(cleanup.base, paste0(sample, "_nucl-blast_db")) |
      startsWith(cleanup.base, paste0(sample, "_last_db"))]
    unlink(c(cleanup.files, red.file, paste0(red.file, ".clstr"), rename.file))

    #Loads in match data
    if (!file.exists(search.file)) stop("Search output was not created for ", sample, ".")
    if (file.size(search.file) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      return(NULL)
    }
    match.data = data.table::fread(search.file, sep = "\t", header = FALSE,
                                   stringsAsFactors = FALSE)
    data.table::setnames(match.data, headers)

    #Matches need to be greater than 12
    filt.data = match.data[match.data$matches > min.match.length,]
    #Percent identitiy must match 50% or greater
    filt.data = filt.data[filt.data$pident >= min.match.percent,]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      return(NULL)
      }

    #Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    # The coverage test is not applied here. It runs after Part C and Part D,
    # once the fragments of a target are joined. See the note there.

    #Reads in contigs
    contigs = all.data

    # Build one oriented candidate per target and contig. Multiple alignments on
    # the same contig describe one copy; separate contigs remain separate copies.
    interval.width = function(starts, ends) {
      intervals = data.frame(start = pmin(starts, ends), end = pmax(starts, ends))
      intervals = intervals[order(intervals$start, intervals$end), , drop = FALSE]
      total = 0
      current.start = intervals$start[1]
      current.end = intervals$end[1]
      if (nrow(intervals) > 1) {
        for (row in 2:nrow(intervals)) {
          if (intervals$start[row] <= current.end + 1) {
            current.end = max(current.end, intervals$end[row])
          } else {
            total = total + current.end - current.start + 1
            current.start = intervals$start[row]
            current.end = intervals$end[row]
          }
        }
      }
      total + current.end - current.start + 1
    }

    fin.loci = Biostrings::DNAStringSet()
    coverage.values = numeric(0)
    target.lengths = numeric(0)

    for (target.name in unique(filt.data$qName)) {
      target.hits = filt.data[filt.data$qName == target.name, ]
      pieces = list()

      for (contig.name in unique(target.hits$tName)) {
        hits = target.hits[target.hits$tName == contig.name, ]
        strand.row = which.max(hits$bitscore)
        same.strand = sign(hits$qEnd[strand.row] - hits$qStart[strand.row]) ==
                      sign(hits$tEnd[strand.row] - hits$tStart[strand.row])
        sequence = contigs[contig.name]
        target.start = min(hits$qStart, hits$qEnd)
        target.end = max(hits$qStart, hits$qEnd)
        contig.start = min(hits$tStart, hits$tEnd)
        contig.end = max(hits$tStart, hits$tEnd)

        if (!same.strand) {
          sequence = Biostrings::reverseComplement(sequence)
          old.start = contig.start
          contig.start = hits$tLen[1] - contig.end + 1
          contig.end = hits$tLen[1] - old.start + 1
        }

        extract.start = max(1, contig.start - (target.start - 1))
        extract.end = min(hits$tLen[1], contig.end + (hits$qLen[1] - target.end))
        pieces[[length(pieces) + 1]] = list(
          sequence = Biostrings::subseq(sequence, start = extract.start,
                                        end = extract.end),
          q.start = target.start,
          q.end = target.end,
          coverage = interval.width(hits$qStart, hits$qEnd),
          target.length = max(hits$qLen)
        )
      }

      piece.order = order(vapply(pieces, `[[`, numeric(1), "q.start"),
                          vapply(pieces, `[[`, numeric(1), "q.end"))
      pieces = pieces[piece.order]
      starts = vapply(pieces, `[[`, numeric(1), "q.start")
      ends = vapply(pieces, `[[`, numeric(1), "q.end")
      distinct.copies = length(pieces) > 1 && any(starts[-1] - ends[-length(ends)] < -30)

      if (distinct.copies) {
        for (piece in pieces) {
          names(piece$sequence) = target.name
          fin.loci = append(fin.loci, piece$sequence)
          coverage.values = c(coverage.values, piece$coverage)
          target.lengths = c(target.lengths, piece$target.length)
        }
      } else {
        joined = character(0)
        for (piece.index in seq_along(pieces)) {
          joined = c(joined, as.character(pieces[[piece.index]]$sequence))
          if (piece.index < length(pieces)) {
            gap = pieces[[piece.index + 1]]$q.start - pieces[[piece.index]]$q.end - 1
            if (gap > 0) joined = c(joined, paste(rep("N", gap), collapse = ""))
          }
        }
        joined.sequence = Biostrings::DNAStringSet(paste(joined, collapse = ""))
        names(joined.sequence) = target.name
        fin.loci = append(fin.loci, joined.sequence)
        coverage.values = c(coverage.values, interval.width(target.hits$qStart,
                                                            target.hits$qEnd))
        target.lengths = c(target.lengths, max(target.hits$qLen))
      }
    }

    keep = .baseWidth(fin.loci) >= min.match.length &
           coverage.values >= ((min.match.coverage / 100) * target.lengths)
    fin.loci = fin.loci[keep]

    if (length(fin.loci) == 0) {
      print(paste0(sample, " had no curated contigs. Skipping"))
      return(NULL)
    }

    # A paralog keeps its own copy. Dropping one here would send its reads to the
    # copy that remains and make false heterozygosity in the variant caller.
    names(fin.loci) = make.unique(names(fin.loci), sep = "_")

    final.loci = as.list(as.character(fin.loci))
    out.file = paste0(output.directory, "/", sample, ".fa")
    temp.file = paste0(out.file, ".tmp-", Sys.getpid())
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      temp.file, nbchar = 1000000, as.string = TRUE
    )
    if (!file.exists(temp.file) || file.size(temp.file) == 0 ||
        file.rename(temp.file, out.file) == FALSE) {
      unlink(temp.file)
      stop("Could not publish curated contigs for ", sample, ".")
    }

    #------------------------------------------------------
    # Per-sample match log
    #------------------------------------------------------
    filt.log = as.data.frame(filt.data)[, c("qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
    filt.log$Sample = sample
    filt.log = filt.log[, c("Sample", "qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
    write.csv(filt.log, file = paste0("logs/sample_logs/", sample, "_curation-matches.csv"), row.names = FALSE)

    print(paste0(sample, " curation complete. ", length(final.loci), " targets kept."))

    return(data.frame(
      Sample          = sample,
      DedupContigs    = length(all.data),
      TargetsMatched  = length(unique(filt.data$qName)),
      CuratedTargets  = length(final.loci),
      MeanPident      = round(mean(filt.data$pident), 2),
      MeanBitscore    = round(mean(filt.data$bitscore), 1),
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    print(paste0(file.names[i], " failed: ", conditionMessage(e)))
    return("failed")
  })
  }, mc.cores = threads)

  # A warning raised in a forked child never reaches the parent, so a failed
  # sample must be counted here.
  fail.count = sum(vapply(results, function(x) identical(x, "failed"), logical(1)))
  if (fail.count != 0) {
    print(paste0(fail.count, " of ", length(file.names),
                 " samples failed. See the messages above."))
  }

  #Writes the cross-sample summary
  summary.df = do.call(rbind, results[vapply(results, is.data.frame, logical(1))])
  if (is.null(summary.df) == FALSE && nrow(summary.df) != 0) {
    write.csv(summary.df, file = "logs/curateTargetContigs_summary.csv", row.names = FALSE)
  }

  return(invisible(NULL))
}#end function
