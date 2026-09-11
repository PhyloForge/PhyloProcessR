#' Collect candidate copies for workflow X5
#'
#' Inventories workflow 4 alignments, candidate tables, saved paralog FASTA
#' files, and optional curated references. It keeps recognizable target and
#' sample names and writes one collected record for each target.
#'
#' @param alignment.directory Directory that contains workflow 4 alignments.
#' @param alignment.format Input alignment format, either `phylip` or `fasta`.
#' @param paralog.directory Directory that contains saved paralog FASTA files.
#' @param primary.directory Directory that contains selected primary FASTA files.
#' @param candidate.directory Directory that contains `*_target-candidates.csv` files.
#' @param output.directory Workflow X5 output directory.
#' @param target.names Optional exact target IDs to process.
#' @param reference.file Optional FASTA file of curated references.
#' @param reference.table Optional TSV file with reference annotations.
#' @param overwrite Replace X5 collection outputs when TRUE.
#'
#' @return A target summary, invisibly.
#' @export
collectParalogCopies = function(
    alignment.directory = "data-analysis/alignments/untrimmed_all-markers",
    alignment.format = c("phylip", "fasta"),
    paralog.directory = "data-analysis/contigs/9_paralog-contigs",
    primary.directory = "data-analysis/contigs/8_annotated-contigs",
    candidate.directory = "logs/sample_logs",
    output.directory = "data-analysis/paralog-analysis",
    target.names = NULL,
    reference.file = NULL,
    reference.table = NULL,
    overwrite = FALSE) {

  alignment.format = match.arg(alignment.format)
  source.paths = c(alignment.directory, paralog.directory, primary.directory,
                   candidate.directory, reference.file, reference.table)
  source.paths = source.paths[!is.na(source.paths) & nzchar(source.paths)]
  .x5CheckOutputPath(output.directory, source.paths)

  extension = if (alignment.format == "phylip") "\\.(phy|phylip)$" else
    "\\.(fa|fas|fasta)$"
  alignment.files = if (dir.exists(alignment.directory)) {
    list.files(alignment.directory, pattern = extension, full.names = TRUE,
               ignore.case = TRUE)
  } else character()
  alignment.targets = sub(extension, "", basename(alignment.files),
                          ignore.case = TRUE)
  if (anyDuplicated(alignment.targets)) {
    stop("More than one base alignment has the same target ID.")
  }
  names(alignment.files) = alignment.targets

  candidate.files = if (dir.exists(candidate.directory)) {
    list.files(candidate.directory, pattern = "_target-candidates\\.csv$",
               full.names = TRUE)
  } else character()
  if (length(candidate.files) == 0 && length(alignment.files) > 0) {
    stop("No workflow 4 candidate tables were found in ", candidate.directory,
         ". Run workflow 4 with candidate logging enabled.")
  }
  candidate.list = lapply(candidate.files, function(file) {
    data = utils::read.csv(file, stringsAsFactors = FALSE,
                           check.names = FALSE)
    missing.columns = setdiff(.x5RequiredCandidateColumns, names(data))
    if (length(missing.columns) > 0) {
      stop("Candidate table ", file, " is missing columns: ",
           paste(missing.columns, collapse = ", "))
    }
    for (column in c("Passes_absolute_filters", "Selected_for_primary",
                     "Saved_as_paralog", "Competitive_with_best")) {
      data[[column]] = .x5CheckLogical(data[[column]], column)
    }
    data$Candidate_file = normalizePath(file)
    data
  })
  candidates = if (length(candidate.list) == 0) {
    data.frame(matrix(nrow = 0, ncol = length(.x5RequiredCandidateColumns)),
               stringsAsFactors = FALSE)
  } else do.call(rbind, candidate.list)
  if (nrow(candidates) > 0) {
    if (any(is.na(candidates$Sample) | candidates$Sample == "") ||
        any(is.na(candidates$Target) | candidates$Target == "") ||
        any(is.na(candidates$Source_contig) |
            candidates$Source_contig == "")) {
      stop("Candidate Sample, Target, and Source_contig values must be nonempty.")
    }
    ranked = !is.na(candidates$Candidate_rank)
    if (any(ranked & (candidates$Candidate_rank < 1 |
                      candidates$Candidate_rank !=
                        as.integer(candidates$Candidate_rank)))) {
      stop("Candidate_rank values must be positive integers when supplied.")
    }
    keys = candidates[candidates$Passes_absolute_filters &
                        !is.na(candidates$Candidate_rank),
                      c("Sample", "Target", "Candidate_rank"), drop = FALSE]
    if (anyDuplicated(keys)) {
      stop("Candidate tables contain duplicate Sample, Target, and Candidate_rank keys.")
    }
    if (any(candidates$Saved_as_paralog &
            (!candidates$Passes_absolute_filters |
             is.na(candidates$Candidate_rank)))) {
      stop("Every saved paralog record must pass the absolute filters and have a rank.")
    }
  }

  saved.targets = sort(unique(candidates$Target[candidates$Saved_as_paralog]))
  all.targets = sort(unique(c(alignment.targets, saved.targets)))
  available.targets = all.targets
  outside.targets = sort(setdiff(unique(candidates$Target), all.targets))
  if (!is.null(target.names)) {
    missing.targets = setdiff(target.names, all.targets)
    if (length(missing.targets) > 0) {
      stop("target.names contains IDs outside the available target union: ",
           paste(missing.targets, collapse = ", "))
    }
    all.targets = sort(unique(target.names))
  }
  if (length(all.targets) == 0) {
    stop("No base alignments or saved paralog targets are available.")
  }

  sequence.files = unique(c(
    if (dir.exists(paralog.directory)) list.files(
      paralog.directory, pattern = "\\.(fa|fas|fasta)$", full.names = TRUE,
      ignore.case = TRUE) else character(),
    if (dir.exists(primary.directory)) list.files(
      primary.directory, pattern = "\\.(fa|fas|fasta)$", full.names = TRUE,
      ignore.case = TRUE) else character()
  ))
  sequence.sets = lapply(sequence.files, .x5ReadSequenceFile)
  names(sequence.sets) = sequence.files
  fasta.index = list()
  for (file in sequence.files) {
    for (header in names(sequence.sets[[file]])) {
      if (!is.null(fasta.index[[header]])) {
        stop("FASTA header occurs in more than one source file: ", header)
      }
      fasta.index[[header]] = list(
        sequence = sequence.sets[[file]][header], file = normalizePath(file)
      )
    }
  }

  references = NULL
  reference.data = NULL
  if (xor(is.null(reference.file), is.null(reference.table))) {
    stop("reference.file and reference.table must be supplied together.")
  }
  if (!is.null(reference.file)) {
    references = .x5ReadSequenceFile(reference.file)
    reference.data = .x5ReadTable(reference.table)
    required = c("Reference_id", "Target", "Role", "Copy_group")
    if (!all(required %in% names(reference.data))) {
      stop("reference.table must contain Reference_id, Target, Role, and Copy_group.")
    }
    if (anyDuplicated(reference.data$Reference_id) ||
        !setequal(reference.data$Reference_id, names(references))) {
      stop("Reference IDs must be unique and match reference.file exactly.")
    }
    if (any(!reference.data$Role %in% c("copy", "outgroup"))) {
      stop("Reference Role must be copy or outgroup.")
    }
    invalid.group = reference.data$Role == "copy" &
      (is.na(reference.data$Copy_group) | reference.data$Copy_group == "")
    if (any(invalid.group)) stop("Copy references require a Copy_group value.")
    invalid.outgroup = reference.data$Role == "outgroup" &
      !is.na(reference.data$Copy_group) & reference.data$Copy_group != ""
    if (any(invalid.outgroup)) stop("Outgroup references cannot have a Copy_group value.")
    if (any(!reference.data$Target %in% available.targets)) {
      stop("reference.table contains a target outside the available target union.")
    }
  }

  input.files = c(alignment.files, candidate.files, sequence.files,
                   reference.file, reference.table)
  settings = list(
    alignment.format = alignment.format,
    target.names = sort(target.names),
    manifest = .x5FileManifest(input.files)
  )
  .x5StageSettings(output.directory, "collect", settings, overwrite)
  .x5InvalidateDownstream(output.directory, "collect", overwrite)
  collected.directory = file.path(output.directory, "1_collected")
  .x5ClearDirectory(collected.directory, overwrite)
  tables.directory = file.path(output.directory, "tables")
  dir.create(tables.directory, recursive = TRUE, showWarnings = FALSE)

  target.rows = list()
  copy.rows = list()
  used.headers = character()
  safe.names = vapply(all.targets, .x5SafeTargetName, character(1))
  if (anyDuplicated(safe.names)) {
    duplicated.names = unique(safe.names[duplicated(safe.names) |
                                         duplicated(safe.names, fromLast = TRUE)])
    stop("Target IDs resolve to the same output file name: ",
         paste(duplicated.names, collapse = ", "))
  }

  for (target.index in seq_along(all.targets)) {
    target = all.targets[target.index]
    locus.id = unname(safe.names[target.index])
    base.file = unname(alignment.files[target])
    has.base = length(base.file) == 1 && !is.na(base.file)
    base = Biostrings::DNAStringSet()
    locus.maps = list()
    additions = Biostrings::DNAStringSet()

    target.candidates = candidates[candidates$Target == target &
                                     candidates$Passes_absolute_filters, ,
                                   drop = FALSE]
    if (has.base) {
      base = .x5ReadAlignment(base.file, alignment.format)
      base = base[order(names(base))]
      if (anyDuplicated(names(base))) {
        stop("Base alignment has duplicate sample labels for target ", target)
      }
      for (sample in names(base)) {
        selected = target.candidates[
          target.candidates$Sample == sample &
            target.candidates$Selected_for_primary, , drop = FALSE]
        if (nrow(selected) != 1) {
          stop("Base sample ", sample, " in target ", target,
               " does not map to one selected workflow 4 candidate.")
        }
        sequence.id = paste0("X5_SEQUENCE_", length(locus.maps) + 1L)
        names(base)[names(base) == sample] = sequence.id
        map = selected[1, , drop = FALSE]
        map$Sequence_id = sequence.id
        map$Source_file = normalizePath(base.file)
        map$Source_header = sample
        map$In_base = TRUE
        map$Added = FALSE
        map$Role = "copy"
        map$Reference_id = ""
        map$Copy_group = ""
        locus.maps[[length(locus.maps) + 1]] = map
      }
    }

    available = target.candidates[
      target.candidates$Saved_as_paralog |
        (!has.base & target.candidates$Selected_for_primary), , drop = FALSE]
    if (nrow(available) > 0) {
      available = available[order(available$Sample, available$Candidate_rank,
                                  available$Source_contig), , drop = FALSE]
      identities = paste(available$Sample, available$Target,
                         available$Candidate_rank, sep = "\r")
      available = available[!duplicated(identities), , drop = FALSE]
    }
    for (row.index in seq_len(nrow(available))) {
      candidate = available[row.index, , drop = FALSE]
      rank = as.integer(candidate$Candidate_rank)
      paralog.header = paste0(target, "_|_", candidate$Sample,
                              "_|_copy", sprintf("%02d", rank))
      primary.header = paste0(target, "_|_", candidate$Sample)
      header = if (candidate$Saved_as_paralog) paralog.header else primary.header
      source = fasta.index[[header]]
      if (is.null(source)) {
        stop("Required FASTA record is missing: ", header)
      }
      used.headers = c(used.headers, header)

      existing = which(vapply(locus.maps, function(map) {
        map$Sample == candidate$Sample &&
          map$Candidate_rank == candidate$Candidate_rank
      }, logical(1)))
      if (length(existing) > 0) {
        existing.id = locus.maps[[existing[1]]]$Sequence_id
        if (!.x5SameSequence(base[existing.id], source$sequence)) {
          stop("Sequence mismatch for candidate ", header,
               " between the base alignment and saved FASTA.")
        }
        next
      }

      sequence.id = paste0("X5_SEQUENCE_", length(locus.maps) + 1L)
      sequence = source$sequence
      names(sequence) = sequence.id
      additions = c(additions, sequence)
      map = candidate
      map$Sequence_id = sequence.id
      map$Source_file = source$file
      map$Source_header = header
      map$In_base = FALSE
      map$Added = TRUE
      map$Role = "copy"
      map$Reference_id = ""
      map$Copy_group = ""
      locus.maps[[length(locus.maps) + 1]] = map
    }

    if (!is.null(reference.data)) {
      target.references = reference.data[reference.data$Target == target, , drop = FALSE]
      for (row.index in seq_len(nrow(target.references))) {
        item = target.references[row.index, , drop = FALSE]
        sequence.id = paste0("X5_SEQUENCE_", length(locus.maps) + 1L)
        sequence = references[item$Reference_id]
        names(sequence) = sequence.id
        additions = c(additions, sequence)
        map = as.data.frame(as.list(rep(NA, length(names(candidates)))),
                            stringsAsFactors = FALSE)
        names(map) = names(candidates)
        map$Sample = ""
        map$Target = target
        map$Candidate_rank = NA_integer_
        map$Sequence_id = sequence.id
        map$Source_file = normalizePath(reference.file)
        map$Source_header = item$Reference_id
        map$In_base = FALSE
        map$Added = TRUE
        map$Role = if (item$Role == "copy") "copy_reference" else
          "outgroup_reference"
        map$Reference_id = item$Reference_id
        map$Copy_group = ifelse(is.na(item$Copy_group), "", item$Copy_group)
        locus.maps[[length(locus.maps) + 1]] = map
      }
    }

    locus.map = do.call(rbind, locus.maps)
    rownames(locus.map) = NULL
    old.ids = locus.map$Sequence_id
    sequence.ids = character(nrow(locus.map))
    biological = locus.map$Role == "copy"
    for (sample in unique(locus.map$Sample[biological])) {
      sample.rows = which(biological & locus.map$Sample == sample)
      order.rows = order(locus.map$Candidate_rank[sample.rows],
                         locus.map$Source_contig[sample.rows], na.last = TRUE)
      sample.rows = sample.rows[order.rows]
      sequence.ids[sample.rows] = if (length(sample.rows) == 1) sample else
        paste0(sample, "_", seq_along(sample.rows))
    }
    reference.rows = which(!biological)
    sequence.ids[reference.rows] = locus.map$Reference_id[reference.rows]
    if (any(is.na(sequence.ids) | sequence.ids == "") ||
        anyDuplicated(sequence.ids)) {
      stop("Sample and reference names are not unique after copy suffixes for target ",
           target, ".")
    }
    names(base) = sequence.ids[match(names(base), old.ids)]
    names(additions) = sequence.ids[match(names(additions), old.ids)]
    locus.map$Sequence_id = sequence.ids
    copy.rows[[target.index]] = locus.map
    record = list(
      Locus_id = locus.id, Target = target, Safe_name = locus.id,
      Base_file = if (has.base) normalizePath(base.file) else "",
      Has_base = has.base, Base_alignment = base, Additions = additions,
      Copy_map = locus.map
    )
    saveRDS(record, file.path(collected.directory, paste0(locus.id, ".rds")))
    biological = locus.map$Role == "copy"
    target.rows[[target.index]] = data.frame(
      Locus_id = locus.id, Target = target,
      Safe_name = locus.id,
      Base_file = if (has.base) normalizePath(base.file) else "",
      Has_base = has.base,
      Initial_samples = length(unique(locus.map$Sample[biological])),
      Initial_copies = sum(biological),
      Additions = sum(locus.map$Added & biological),
      Newly_constructed = !has.base,
      Run_scope = if (is.null(target.names)) "full" else "subset",
      stringsAsFactors = FALSE
    )
  }

  target.map = do.call(rbind, target.rows)
  copy.map = do.call(rbind, copy.rows)
  copy.map = copy.map[, c("Sequence_id", "Target", "Sample",
                          "Candidate_rank", "Source_contig", "Source_file",
                          "Source_header", "In_base", "Added", "Role",
                          "Reference_id", "Copy_group",
                          setdiff(names(copy.map), c(
                            "Sequence_id", "Target", "Sample", "Candidate_rank",
                            "Source_contig", "Source_file", "Source_header",
                            "In_base", "Added", "Role", "Reference_id",
                            "Copy_group"))), drop = FALSE]
  utils::write.table(target.map, file.path(tables.directory, "target-map.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(copy.map, file.path(tables.directory, "copy-map.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  unreferenced = setdiff(names(fasta.index), used.headers)
  utils::write.table(data.frame(Source_header = unreferenced),
                     file.path(tables.directory, "unreferenced-sources.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(data.frame(Target = outside.targets),
                     file.path(tables.directory, "candidate-targets-outside-union.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  saveRDS(list(stage = "collect", complete = TRUE, target.map = target.map),
          file.path(collected.directory, ".complete.rds"))
  message(nrow(target.map), " targets collected. Tables: ", tables.directory)
  invisible(target.map)
}
