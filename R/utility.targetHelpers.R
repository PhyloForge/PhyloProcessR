# Internal helpers that assign contigs to targets and compare copies of one
# target. Used by curateTargetContigs and annotateTargets.

# Removes the weaker target where two targets match the same region of one
# contig. Two spans share a region when they overlap by more than half of the
# shorter span. The span with the higher score is kept. When the scores are
# equal, the target name that sorts first is kept.
.keepBestTargetPerRegion = function(spans) {
  spans = spans[order(spans$tName, -spans$score, spans$qName)]
  keep = rep(TRUE, nrow(spans))
  rows.by.contig = split(seq_len(nrow(spans)), spans$tName)

  for (rows in rows.by.contig) {
    if (length(rows) < 2) next
    # Compare each span with the stronger spans on the same contig
    for (i in 2:length(rows)) {
      for (j in 1:(i - 1)) {
        if (keep[rows[j]] == FALSE) next
        start = max(spans$target.start[rows[i]], spans$target.start[rows[j]])
        end = min(spans$target.end[rows[i]], spans$target.end[rows[j]])
        overlap = end - start + 1
        shorter = min(spans$target.end[rows[i]] - spans$target.start[rows[i]] + 1,
                      spans$target.end[rows[j]] - spans$target.start[rows[j]] + 1)
        if (overlap > shorter / 2) {
          keep[rows[i]] = FALSE
          break
        }
      }
    }
  }
  spans[keep]
}

# Drops the pieces of one target that are the same copy as a longer piece:
# their overlap is at least min.identity percent identical. Pieces are checked
# from the longest to the shortest, and a piece is kept only when it differs
# from every piece kept before it.
.dropSameCopyPieces = function(pieces, min.identity, min.columns) {
  sequences = lapply(pieces, function(piece) {
    Biostrings::DNAString(paste0(as.character(piece$left.flank),
                                 as.character(piece$sequence),
                                 as.character(piece$right.flank)))
  })
  longest.first = order(-vapply(sequences, length, numeric(1)))
  kept = integer(0)

  for (index in longest.first) {
    same.copy = FALSE
    for (kept.index in kept) {
      identity = .overlapIdentity(sequences[[kept.index]], sequences[[index]],
                                  min.columns = min.columns)
      if (!is.na(identity) && identity >= min.identity) {
        same.copy = TRUE
        break
      }
    }
    if (same.copy == FALSE) kept = c(kept, index)
  }
  pieces[sort(kept)]
}

# Percent identity of two contigs where they overlap. The contigs are aligned
# end to end with no penalty for the overhangs, in both orientations. Only
# columns where both contigs have A, C, G or T are counted. Returns NA when no
# orientation gives at least min.columns such columns.
.overlapIdentity = function(first, second, min.columns = 50) {
  # The pairwise alignment functions moved from Biostrings to pwalign in
  # Bioconductor 3.19. Use pwalign when it is installed; older Biostrings
  # versions still have the functions.
  package = if (requireNamespace("pwalign", quietly = TRUE)) "pwalign" else "Biostrings"
  pairwise.alignment = getExportedValue(package, "pairwiseAlignment")
  aligned.pattern = getExportedValue(package, "alignedPattern")
  aligned.subject = getExportedValue(package, "alignedSubject")
  substitution.matrix = getExportedValue(package, "nucleotideSubstitutionMatrix")

  scoring = substitution.matrix(match = 1, mismatch = -3, baseOnly = FALSE)
  best.identity = NA_real_
  for (subject in list(second, Biostrings::reverseComplement(second))) {
    alignment = pairwise.alignment(first, subject, type = "overlap",
                                   substitutionMatrix = scoring,
                                   gapOpening = 5, gapExtension = 2)
    first.bases = strsplit(as.character(aligned.pattern(alignment)), "")[[1]]
    second.bases = strsplit(as.character(aligned.subject(alignment)), "")[[1]]
    both = first.bases %in% c("A", "C", "G", "T") & second.bases %in% c("A", "C", "G", "T")
    if (sum(both) < min.columns) next
    identity = 100 * sum(first.bases[both] == second.bases[both]) / sum(both)
    if (is.na(best.identity) || identity > best.identity) best.identity = identity
  }
  best.identity
}
