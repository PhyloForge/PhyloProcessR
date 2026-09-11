#' Expand alignments with candidate paralog copies
#'
#' Preserves existing workflow 4 alignments and adds missing copies with MAFFT.
#' It creates a new alignment when a saved-copy target has no base alignment.
#'
#' @param collected.directory Directory made by [collectParalogCopies()].
#' @param target.file Capture-target FASTA file used to guide new alignments.
#' @param output.directory Workflow X5 output directory.
#' @param alignment.algorithm MAFFT algorithm for new alignments.
#' @param threads Number of MAFFT threads.
#' @param mafft.path Directory that contains MAFFT, or NULL for PATH.
#' @param quiet Suppress MAFFT console output when TRUE.
#' @param overwrite Replace alignment expansion outputs when TRUE.
#'
#' @return An alignment expansion summary, invisibly.
#' @export
alignParalogCopies = function(
    collected.directory = "data-analysis/paralog-analysis/1_collected",
    target.file = NULL,
    output.directory = "data-analysis/paralog-analysis",
    alignment.algorithm = c("localpair", "globalpair"),
    threads = 8,
    mafft.path = NULL,
    quiet = TRUE,
    overwrite = FALSE) {

  alignment.algorithm = match.arg(alignment.algorithm)
  .validateResources(threads = threads, memory = 1)
  records = list.files(collected.directory, pattern = "^[^.].*\\.rds$",
                       full.names = TRUE)
  if (length(records) == 0) stop("No collected X5 target records were found.")
  if (is.null(target.file) || !file.exists(target.file)) {
    if (any(vapply(records, function(file) !readRDS(file)$Has_base, logical(1)))) {
      stop("target.file is required when a target has no base alignment.")
    }
  }
  manifest.paths = c(records, target.file)
  settings = list(
    alignment.algorithm = alignment.algorithm, threads = as.integer(threads),
    mafft.path = mafft.path, manifest = .x5FileManifest(manifest.paths)
  )
  .x5StageSettings(output.directory, "align", settings, overwrite)
  .x5InvalidateDownstream(output.directory, "align", overwrite)
  expanded.directory = file.path(output.directory, "2_expanded")
  work.root = file.path(output.directory, "logs", "alignment")
  dir.create(file.path(output.directory, "tables"), recursive = TRUE,
             showWarnings = FALSE)
  .x5ClearDirectory(expanded.directory, overwrite)
  .x5ClearDirectory(work.root, overwrite)

  guides = if (!is.null(target.file) && file.exists(target.file)) {
    .x5ReadSequenceFile(target.file)
  } else NULL
  executable = NULL
  rows = list()

  run.mafft = function(command, output, log, task) {
    unlink(c(output, log))
    .runCommand(paste(command, ">", shQuote(output)), quiet = quiet,
                task = task, keep.stdout = TRUE, stderr.log = log)
    .x5ReadAlignment(output, "fasta")
  }

  for (index in seq_along(records)) {
    record = readRDS(records[index])
    result.file = file.path(expanded.directory,
                            paste0(record$Locus_id, ".rds"))
    fasta.file = file.path(expanded.directory,
                           paste0(record$Locus_id, ".fa"))
    if (!overwrite && file.exists(result.file) && file.exists(fasta.file)) {
      result = readRDS(result.file)
      .x5ValidateAlignment(.x5ReadAlignment(fasta.file, "fasta"),
                           record$Locus_id, result$Sequence_ids)
      rows[[index]] = result$summary
      next
    }

    base = record$Base_alignment
    additions = record$Additions
    path = if (record$Has_base && length(additions) == 0) "unchanged_base" else
      if (record$Has_base) "added_copies" else "new_alignment"
    orientation.changed = character()
    work.directory = file.path(work.root, record$Locus_id)
    dir.create(work.directory, recursive = TRUE, showWarnings = FALSE)
    log.file = file.path(work.directory, "mafft.log")

    if (path == "unchanged_base") {
      alignment = base
    } else {
      if (is.null(executable)) executable = .toolCommand("mafft", mafft.path)
      input.file = file.path(work.directory, "base.fa")
      additions.file = file.path(work.directory, "additions.fa")
      output.file = file.path(work.directory, "alignment.fa")
      if (path == "added_copies") {
        .x5WriteFasta(base, input.file)
        .x5WriteFasta(additions, additions.file)
        command = paste(
          executable, "--add", shQuote(additions.file), "--maxiterate 0",
          "--nuc --adjustdirection --thread", as.integer(threads),
          shQuote(input.file)
        )
        alignment = run.mafft(command, output.file, log.file,
                              paste(record$Target, "MAFFT add"))
      } else {
        if (!record$Target %in% names(guides)) {
          stop("target.file has no exact guide for target ", record$Target)
        }
        guide = guides[record$Target]
        names(guide) = "X5_GUIDE"
        input = c(additions, guide)
        .x5WriteFasta(input, input.file)
        command = paste(
          executable, paste0("--", alignment.algorithm),
          "--maxiterate 1000 --nuc --adjustdirection --quiet --op 3",
          "--ep 0.123 --thread", as.integer(threads), shQuote(input.file)
        )
        alignment = run.mafft(command, output.file, log.file,
                              paste(record$Target, "MAFFT alignment"))
        names(alignment) = sub("^_R_", "", names(alignment))
        alignment = alignment[names(alignment) != "X5_GUIDE"]
      }
    }

    original.names = names(alignment)
    normalized.names = sub("^_R_", "", original.names)
    orientation.changed = normalized.names[normalized.names != original.names]
    names(alignment) = normalized.names
    expected = c(names(base), names(additions))
    alignment = .x5ValidateAlignment(alignment, record$Locus_id, expected)

    source.sequences = c(base, additions)
    for (sequence.id in expected) {
      if (!.x5SameSequence(alignment[sequence.id], source.sequences[sequence.id])) {
        stop("MAFFT changed the ungapped sequence content for ", sequence.id)
      }
    }
    if (record$Has_base && path == "added_copies") {
      base.output = alignment[names(base)]
      base.matrix = do.call(rbind, strsplit(as.character(base.output), ""))
      keep = colSums(base.matrix != "-") > 0
      base.matrix = base.matrix[, keep, drop = FALSE]
      original.matrix = do.call(rbind, strsplit(as.character(base), ""))
      if (!identical(unname(base.matrix), unname(original.matrix))) {
        stop("MAFFT changed base-alignment homology for target ", record$Target)
      }
    }

    .x5WriteFasta(alignment, fasta.file)
    summary = data.frame(
      Locus_id = record$Locus_id, Target = record$Target,
      Expansion_path = path, Sequence_count = length(alignment),
      Alignment_length = Biostrings::width(alignment)[1],
      Orientation_changes = length(orientation.changed),
      Alignment_file = normalizePath(fasta.file), stringsAsFactors = FALSE
    )
    result = list(
      stage = "align", complete = TRUE, Target = record$Target,
      Locus_id = record$Locus_id, Sequence_ids = names(alignment),
      orientation.changed = orientation.changed, summary = summary
    )
    saveRDS(result, result.file)
    rows[[index]] = summary
  }
  summary = do.call(rbind, rows)
  utils::write.table(summary,
                     file.path(output.directory, "tables", "alignment-expansion.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  message(nrow(summary), " targets expanded. Alignments: ", expanded.directory)
  invisible(summary)
}
