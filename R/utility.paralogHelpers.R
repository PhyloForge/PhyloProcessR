# Internal helpers for workflow X5 paralog analysis.

.x5RequiredCandidateColumns = c(
  "Sample", "Target", "Contig", "Source_contig", "Candidate_rank",
  "Passes_absolute_filters", "Selected_for_primary", "Saved_as_paralog",
  "Competitive_with_best"
)

.x5CheckLogical = function(values, column) {
  text = toupper(as.character(values))
  valid = !is.na(text) & text %in% c("TRUE", "FALSE")
  if (!all(valid)) {
    stop(column, " must contain only TRUE or FALSE values.")
  }
  text == "TRUE"
}

.x5ReadAlignment = function(file, format = c("phylip", "fasta")) {
  format = match.arg(format)
  if (!file.exists(file) || file.info(file)$size == 0) {
    stop("Alignment file is missing or empty: ", file)
  }
  if (format == "phylip") {
    alignment = Biostrings::DNAStringSet(
      Biostrings::readDNAMultipleAlignment(file, format = "phylip")
    )
  } else {
    alignment = Biostrings::readDNAStringSet(file)
  }
  .x5ValidateAlignment(alignment, basename(file))
}

.x5ValidateAlignment = function(alignment, label = "alignment",
                                 expected.ids = NULL) {
  if (!inherits(alignment, "DNAStringSet") || length(alignment) == 0) {
    stop(label, " is empty or is not a DNAStringSet.")
  }
  if (is.null(names(alignment)) || any(names(alignment) == "") ||
      anyDuplicated(names(alignment))) {
    stop(label, " must have unique, nonempty sequence names.")
  }
  widths = Biostrings::width(alignment)
  if (any(widths == 0) || length(unique(widths)) != 1) {
    stop(label, " must be a nonempty rectangular alignment.")
  }
  if (!is.null(expected.ids) &&
      !setequal(names(alignment), expected.ids)) {
    stop(label, " does not contain the expected sequence IDs.")
  }
  alignment
}

.x5WriteFasta = function(alignment, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary = tempfile("x5-fasta-", tmpdir = dirname(file))
  on.exit(unlink(temporary), add = TRUE)
  writeFasta(
    sequences = as.list(as.character(alignment)),
    names = names(alignment), file.out = temporary,
    nbchar = 1000000, as.string = TRUE
  )
  if (!file.rename(temporary, file)) {
    stop("Could not publish FASTA file: ", file)
  }
  invisible(file)
}

.x5WritePhylip = function(alignment, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  temporary = tempfile("x5-phylip-", tmpdir = dirname(file))
  on.exit(unlink(temporary), add = TRUE)
  matrix.alignment = do.call(rbind, strsplit(as.character(alignment), ""))
  rownames(matrix.alignment) = names(alignment)
  writePhylip(matrix.alignment, file = temporary, interleave = FALSE)
  checked = .x5ReadAlignment(temporary, format = "phylip")
  if (!identical(names(checked), names(alignment))) {
    stop("Sample names did not pass a relaxed PHYLIP round trip for ", file)
  }
  if (!file.rename(temporary, file)) {
    stop("Could not publish PHYLIP file: ", file)
  }
  invisible(file)
}

.x5Ungapped = function(sequence) {
  unname(toupper(gsub("[-?.]", "", as.character(sequence))))
}

.x5SameSequence = function(first, second) {
  first = .x5Ungapped(first)
  second = .x5Ungapped(second)
  reverse = as.character(Biostrings::reverseComplement(
    Biostrings::DNAString(second)
  ))
  identical(first, second) || identical(first, reverse)
}

.x5ReadSequenceFile = function(file) {
  if (!file.exists(file) || file.info(file)$size == 0) {
    return(Biostrings::DNAStringSet())
  }
  sequences = Biostrings::readDNAStringSet(file)
  if (is.null(names(sequences)) || anyDuplicated(names(sequences))) {
    stop("FASTA headers must be unique in ", file)
  }
  sequences
}

.x5SafeTargetName = function(target) {
  safe = gsub("[^A-Za-z0-9_.-]", "_", target)
  safe = sub("^\\.+", "", safe)
  if (safe == "") safe = "target"
  safe
}

.x5CheckOutputPath = function(output.directory, source.paths) {
  output = normalizePath(output.directory, mustWork = FALSE)
  sources = normalizePath(source.paths, mustWork = FALSE)
  inside = sources == output | startsWith(sources, paste0(output, "/")) |
    startsWith(output, paste0(sources, "/"))
  if (any(inside)) {
    stop("The X5 output directory must be separate from all input paths.")
  }
  invisible(output)
}

.x5FileManifest = function(paths) {
  paths = sort(unique(paths[file.exists(paths)]))
  data.frame(
    Path = normalizePath(paths),
    MD5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

.x5StageSettings = function(output.directory, stage, settings,
                            overwrite = FALSE) {
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  file = file.path(output.directory, paste0(".", stage, "-settings.rds"))
  if (file.exists(file)) {
    previous = readRDS(file)
    if (!identical(previous, settings) && !overwrite) {
      stop("Inputs or settings changed for X5 stage ", stage,
           ". Use overwrite = TRUE or a new output directory.")
    }
  }
  saveRDS(settings, file)
  invisible(file)
}

.x5ClearDirectory = function(directory, overwrite = FALSE) {
  if (overwrite && dir.exists(directory)) unlink(directory, recursive = TRUE)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  invisible(directory)
}

.x5InvalidateDownstream = function(output.directory, stage,
                                    overwrite = FALSE) {
  if (!overwrite) return(invisible(NULL))
  paths = switch(
    stage,
    collect = c("2_expanded", "3_tree-alignments", "4_trees", "5_separated",
                "trimmed_all-markers", "retained-fasta", "excluded",
                "provisional-groups", "trimmed_genes", "trimmed_unlinked_genes",
                "trimmed_all-unlinked",
                ".align-settings.rds", ".trim-settings.rds",
                ".trees-settings.rds", ".separate-settings.rds"),
    align = c("3_tree-alignments", "4_trees", "5_separated",
              "trimmed_all-markers", "retained-fasta", "excluded",
              "provisional-groups", "trimmed_genes", "trimmed_unlinked_genes",
              "trimmed_all-unlinked",
              ".trim-settings.rds", ".trees-settings.rds",
              ".separate-settings.rds"),
    trim = c("4_trees", "5_separated", "trimmed_all-markers",
             "retained-fasta", "excluded", "provisional-groups",
             "trimmed_genes", "trimmed_unlinked_genes", "trimmed_all-unlinked",
             ".trees-settings.rds", ".separate-settings.rds"),
    trees = c("5_separated", "trimmed_all-markers", "retained-fasta",
              file.path("excluded", "untrimmed"), "provisional-groups",
              "trimmed_genes", "trimmed_unlinked_genes", "trimmed_all-unlinked",
              ".separate-settings.rds"),
    separate = c("trimmed_genes", "trimmed_unlinked_genes",
                 "trimmed_all-unlinked"),
    character()
  )
  paths = file.path(output.directory, paths)
  existing = paths[file.exists(paths) | dir.exists(paths)]
  if (length(existing) > 0) unlink(existing, recursive = TRUE)
  invisible(NULL)
}

.x5PresentMatrix = function(alignment) {
  matrix = do.call(rbind, strsplit(toupper(as.character(alignment)), ""))
  rownames(matrix) = names(alignment)
  present = matrix %in% c("A", "C", "G", "T", "R", "Y", "S", "W", "K", "M",
                          "B", "D", "H", "V")
  dim(present) = dim(matrix)
  dimnames(present) = dimnames(matrix)
  present
}

.x5RunTrimal = function(alignment, work.directory, trimAl.path = NULL,
                        quiet = TRUE) {
  dir.create(work.directory, recursive = TRUE, showWarnings = FALSE)
  input = file.path(work.directory, "input.fa")
  output = file.path(work.directory, "output.fa")
  log = file.path(work.directory, "trimal.log")
  unlink(c(output, log))
  .x5WriteFasta(alignment, input)
  executable = .toolCommand("trimal", trimAl.path)
  command = paste(
    executable, "-in", shQuote(input), "-out", shQuote(output),
    "-automated1"
  )
  .runCommand(command, quiet = quiet, task = "TrimAl", stderr.log = log)
  trimmed = .x5ReadAlignment(output, "fasta")
  .x5ValidateAlignment(trimmed, "TrimAl output", names(alignment))
}

.x5TrimCopies = function(alignment, copy.map, run.TrimAl = TRUE,
                          trimAl.path = NULL, work.directory = tempdir(),
                          min.external.percent = 50,
                          min.column.gap.percent = 50,
                          min.coverage.percent = 35,
                          min.coverage.bp = 60,
                          min.alignment.length = 100,
                          min.taxa.alignment = 4,
                          max.alignment.gap.percent = 50,
                          quiet = TRUE) {
  alignment = .x5ValidateAlignment(alignment)
  copy.map = copy.map[match(names(alignment), copy.map$Sequence_id), , drop = FALSE]
  if (any(is.na(copy.map$Sequence_id))) {
    stop("The copy map does not contain all alignment sequence IDs.")
  }
  removed = data.frame(
    Sequence_id = character(), Stage = character(), Reason = character(),
    Present_bp = integer(), Coverage_percent = numeric(),
    stringsAsFactors = FALSE
  )

  if (run.TrimAl) {
    alignment = .x5RunTrimal(alignment, work.directory, trimAl.path, quiet)
    alignment = alignment[names(alignment)]
    copy.map = copy.map[match(names(alignment), copy.map$Sequence_id), , drop = FALSE]
  }

  biological = copy.map$Role == "copy"
  sample.names = sort(unique(copy.map$Sample[biological]))
  present = .x5PresentMatrix(alignment)
  sample.present = vapply(sample.names, function(sample) {
    colSums(present[biological & copy.map$Sample == sample, , drop = FALSE]) > 0
  }, logical(ncol(present)))
  if (length(sample.names) == 1) {
    sample.present = matrix(sample.present, ncol = 1)
  }
  occupancy = rowMeans(sample.present) * 100
  edge.columns = which(occupancy >= min.external.percent)
  if (length(edge.columns) == 0) {
    return(list(pass = FALSE, reason = "failed_qc", alignment = NULL,
                removed = removed, sample.count = length(sample.names),
                copy.count = sum(biological), length = 0,
                gap.percent = 100))
  }
  keep = seq.int(min(edge.columns), max(edge.columns))
  alignment = Biostrings::subseq(alignment, start = min(keep), end = max(keep))
  present = .x5PresentMatrix(alignment)

  sample.present = vapply(sample.names, function(sample) {
    colSums(present[biological & copy.map$Sample == sample, , drop = FALSE]) > 0
  }, logical(ncol(present)))
  if (length(sample.names) == 1) {
    sample.present = matrix(sample.present, ncol = 1)
  }
  missing.percent = (1 - rowMeans(sample.present)) * 100
  keep.columns = missing.percent < min.column.gap.percent
  if (!any(keep.columns)) {
    return(list(pass = FALSE, reason = "failed_qc", alignment = NULL,
                removed = removed, sample.count = length(sample.names),
                copy.count = sum(biological), length = 0,
                gap.percent = 100))
  }
  sequence.matrix = do.call(rbind, strsplit(as.character(alignment), ""))
  sequence.matrix = sequence.matrix[, keep.columns, drop = FALSE]
  alignment = Biostrings::DNAStringSet(apply(sequence.matrix, 1, paste0,
                                              collapse = ""))
  names(alignment) = rownames(sequence.matrix)

  present = .x5PresentMatrix(alignment)
  present.bp = rowSums(present)
  coverage = present.bp / ncol(present) * 100
  keep.rows = present.bp >= min.coverage.bp &
    coverage >= min.coverage.percent
  if (any(!keep.rows)) {
    removed = rbind(removed, data.frame(
      Sequence_id = names(alignment)[!keep.rows],
      Stage = "copy_coverage", Reason = "low_copy_coverage",
      Present_bp = present.bp[!keep.rows],
      Coverage_percent = coverage[!keep.rows], stringsAsFactors = FALSE
    ))
  }
  alignment = alignment[keep.rows]
  copy.map = copy.map[keep.rows, , drop = FALSE]
  if (length(alignment) == 0 || !any(copy.map$Role == "copy")) {
    return(list(pass = FALSE, reason = "failed_qc", alignment = NULL,
                removed = removed, sample.count = 0, copy.count = 0,
                length = 0, gap.percent = 100))
  }

  present = .x5PresentMatrix(alignment)
  biological = copy.map$Role == "copy"
  keep.columns = colSums(present[biological, , drop = FALSE]) > 0
  if (!any(keep.columns)) {
    return(list(pass = FALSE, reason = "failed_qc", alignment = NULL,
                removed = removed, sample.count = 0, copy.count = 0,
                length = 0, gap.percent = 100))
  }
  sequence.matrix = do.call(rbind, strsplit(as.character(alignment), ""))
  sequence.matrix = sequence.matrix[, keep.columns, drop = FALSE]
  alignment = Biostrings::DNAStringSet(apply(sequence.matrix, 1, paste0,
                                              collapse = ""))
  names(alignment) = rownames(sequence.matrix)
  present = .x5PresentMatrix(alignment)
  biological = copy.map$Role == "copy"
  sample.count = length(unique(copy.map$Sample[biological]))
  copy.count = sum(biological)
  alignment.length = Biostrings::width(alignment)[1]
  gap.percent = (1 - mean(present[biological, , drop = FALSE])) * 100
  pass = sample.count >= min.taxa.alignment &&
    alignment.length >= min.alignment.length &&
    gap.percent <= max.alignment.gap.percent
  list(
    pass = pass, reason = if (pass) "passed" else "failed_qc",
    alignment = alignment, removed = removed,
    sample.count = sample.count, copy.count = copy.count,
    length = alignment.length, gap.percent = gap.percent
  )
}

.x5VariableSite = function(alignment, copy.map) {
  biological = copy.map$Role[match(names(alignment), copy.map$Sequence_id)] == "copy"
  matrix = do.call(rbind, strsplit(toupper(as.character(alignment[biological])), ""))
  any(vapply(seq_len(ncol(matrix)), function(column) {
    length(unique(matrix[, column][matrix[, column] %in% c("A", "C", "G", "T")])) > 1
  }, logical(1)))
}

.x5ValidateTree = function(tree, expected.ids, label = "tree") {
  if (!inherits(tree, "phylo") ||
      !setequal(tree$tip.label, expected.ids)) {
    stop(label, " has an invalid tip set.")
  }
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length < 0)) {
    stop(label, " has invalid branch lengths.")
  }
  support = suppressWarnings(as.numeric(tree$node.label))
  supplied = !is.na(support)
  if (any(supplied & (support < 0 | support > 100))) {
    stop(label, " has support values outside 0 to 100.")
  }
  tree
}

.x5CanonicalSide = function(first, second) {
  first.key = paste(sort(first), collapse = "|")
  second.key = paste(sort(second), collapse = "|")
  if (first.key <= second.key) list(a = sort(first), b = sort(second)) else
    list(a = sort(second), b = sort(first))
}

.x5TreeSplits = function(tree, copy.map, min.branch.support = 95,
                          min.shared.samples = 3,
                          min.shared.sample.fraction = 0.5,
                          min.split.branch.length = 0.05,
                          min.split.branch.ratio = 10,
                          min.taxa.alignment = 4) {
  tip.count = length(tree$tip.label)
  all.tips = tree$tip.label
  biological.samples = copy.map$Sample[copy.map$Role == "copy"]
  has.multiple.copies = anyDuplicated(biological.samples) > 0
  descendants = function(node) {
    children = tree$edge[tree$edge[, 1] == node, 2]
    tips = children[children <= tip.count]
    internal = children[children > tip.count]
    for (child in internal) tips = c(tips, descendants(child))
    tips
  }
  rows = list()
  internal.lengths = tree$edge.length[tree$edge[, 2] > tip.count]
  median.length = stats::median(internal.lengths[is.finite(internal.lengths) &
                                                   internal.lengths > 0])
  if (!is.finite(median.length)) median.length = NA_real_
  for (edge.index in which(tree$edge[, 2] > tip.count)) {
    node = tree$edge[edge.index, 2]
    first = all.tips[descendants(node)]
    second = setdiff(all.tips, first)
    if (length(first) < 2 || length(second) < 2) next
    sides = .x5CanonicalSide(first, second)
    label.index = node - tip.count
    support = if (is.null(tree$node.label) ||
                  label.index > length(tree$node.label)) NA_real_ else
      suppressWarnings(as.numeric(tree$node.label[label.index]))
    length.value = tree$edge.length[edge.index]
    rows[[length(rows) + 1]] = data.frame(
      Split_id = paste(sides$a, collapse = "|"),
      Side_a = I(list(sides$a)), Side_b = I(list(sides$b)),
      Support = support, Length = length.value,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(data.frame())
  splits = do.call(rbind, rows)
  keys = unique(splits$Split_id)
  results = vector("list", length(keys))
  for (index in seq_along(keys)) {
    matches = which(splits$Split_id == keys[index])
    side.a = splits$Side_a[[matches[1]]]
    side.b = splits$Side_b[[matches[1]]]
    edge.length = sum(splits$Length[matches], na.rm = TRUE)
    supports = splits$Support[matches]
    support = if (all(is.na(supports))) NA_real_ else max(supports, na.rm = TRUE)
    map.a = copy.map[match(side.a, copy.map$Sequence_id), , drop = FALSE]
    map.b = copy.map[match(side.b, copy.map$Sequence_id), , drop = FALSE]
    samples.a = map.a$Sample[map.a$Role == "copy"]
    samples.b = map.b$Sample[map.b$Role == "copy"]
    unique.a = unique(samples.a)
    unique.b = unique(samples.b)
    shared = intersect(unique.a, unique.b)
    denominator = min(length(unique.a), length(unique.b))
    shared.fraction = if (denominator == 0) 0 else length(shared) / denominator
    repeated.a = anyDuplicated(samples.a) > 0
    repeated.b = anyDuplicated(samples.b) > 0
    reference.conflict = .x5ReferenceConflict(map.a, map.b)
    length.ratio = if (is.na(median.length)) NA_real_ else edge.length / median.length
    deep = has.multiple.copies &&
      !is.na(support) && support >= min.branch.support &&
      length(unique.a) >= min.taxa.alignment &&
      length(unique.b) >= min.taxa.alignment &&
      edge.length >= min.split.branch.length &&
      !is.na(length.ratio) && length.ratio >= min.split.branch.ratio
    eligible = !is.na(support) && support >= min.branch.support &&
      edge.length > 0 && length(unique.a) >= min.taxa.alignment &&
      length(unique.b) >= min.taxa.alignment && !repeated.a && !repeated.b &&
      length(shared) >= min.shared.samples &&
      shared.fraction >= min.shared.sample.fraction && !reference.conflict
    failures = character()
    if (is.na(support) || support < min.branch.support) failures = c(failures, "support")
    if (!is.finite(edge.length) || edge.length <= 0) failures = c(failures, "branch_length")
    if (length(unique.a) < min.taxa.alignment || length(unique.b) < min.taxa.alignment) {
      failures = c(failures, "sample_count")
    }
    if (repeated.a || repeated.b) failures = c(failures, "within_group_duplicates")
    if (length(shared) < min.shared.samples) failures = c(failures, "shared_sample_count")
    if (shared.fraction < min.shared.sample.fraction) {
      failures = c(failures, "shared_sample_fraction")
    }
    if (reference.conflict) failures = c(failures, "reference_conflict")
    results[[index]] = data.frame(
      Split_id = keys[index], Support = support, Length = edge.length,
      Length_ratio = length.ratio, Side_a_copy_count = sum(map.a$Role == "copy"),
      Side_b_copy_count = sum(map.b$Role == "copy"),
      Side_a_sample_count = length(unique.a),
      Side_b_sample_count = length(unique.b),
      Side_a_repeated_samples = sum(duplicated(samples.a)),
      Side_b_repeated_samples = sum(duplicated(samples.b)),
      Shared_sample_count = length(shared), Shared_sample_fraction = shared.fraction,
      Shared_samples = paste(sort(shared), collapse = ";"),
      Deep_split_review = deep, Eligible = eligible,
      Failure_reason = paste(failures, collapse = ";"),
      Side_a = I(list(side.a)), Side_b = I(list(side.b)),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, results)
}

.x5ReferenceConflict = function(first, second) {
  first = first[first$Role == "copy_reference", , drop = FALSE]
  second = second[second$Role == "copy_reference", , drop = FALSE]
  labels.first = unique(first$Copy_group[!is.na(first$Copy_group) &
                                          nzchar(first$Copy_group)])
  labels.second = unique(second$Copy_group[!is.na(second$Copy_group) &
                                            nzchar(second$Copy_group)])
  length(labels.first) > 1 || length(labels.second) > 1 ||
    length(intersect(labels.first, labels.second)) > 0
}

.x5ReadTable = function(file) {
  if (is.null(file)) return(NULL)
  as.data.frame(data.table::fread(file, sep = "\t", header = TRUE,
                                  data.table = FALSE))
}
