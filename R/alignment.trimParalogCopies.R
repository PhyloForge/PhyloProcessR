#' Trim expanded paralog-copy alignments
#'
#' Applies copy-aware edge, column, and sequence-coverage filters. Biological
#' samples contribute once to occupancy even when they have multiple copies.
#'
#' @param expanded.directory Directory made by [alignParalogCopies()].
#' @param collected.directory Directory made by [collectParalogCopies()].
#' @param output.directory Workflow X5 output directory.
#' @param run.TrimAl Run TrimAl with `-automated1` before native filters.
#' @param trimAl.path Directory that contains TrimAl, or NULL for PATH.
#' @param min.external.percent Minimum unique-sample occupancy at alignment edges.
#' @param min.column.gap.percent Missing-sample percentage that removes a column.
#' @param min.coverage.percent Minimum present-base percentage for each copy.
#' @param min.coverage.bp Minimum present-base count for each copy.
#' @param min.alignment.length Minimum retained alignment width.
#' @param min.taxa.alignment Minimum retained biological sample count.
#' @param max.alignment.gap.percent Maximum missing percentage in biological copies.
#' @param quiet Suppress external-tool console output when TRUE.
#' @param overwrite Replace X5 trimming outputs when TRUE.
#'
#' @return A trimming summary, invisibly.
#' @export
trimParalogCopies = function(
    expanded.directory = "data-analysis/paralog-analysis/2_expanded",
    collected.directory = "data-analysis/paralog-analysis/1_collected",
    output.directory = "data-analysis/paralog-analysis",
    run.TrimAl = TRUE,
    trimAl.path = NULL,
    min.external.percent = 50,
    min.column.gap.percent = 50,
    min.coverage.percent = 35,
    min.coverage.bp = 60,
    min.alignment.length = 100,
    min.taxa.alignment = 4,
    max.alignment.gap.percent = 50,
    quiet = TRUE,
    overwrite = FALSE) {

  percentages = c(min.external.percent, min.column.gap.percent,
                  min.coverage.percent, max.alignment.gap.percent)
  if (any(!is.finite(percentages)) || any(percentages < 0) ||
      any(percentages > 100)) stop("Trimming percentages must be from 0 to 100.")
  counts = c(min.coverage.bp, min.alignment.length, min.taxa.alignment)
  if (any(!is.finite(counts)) || any(counts < 0)) {
    stop("Trimming count thresholds must be nonnegative.")
  }
  alignment.files = list.files(expanded.directory,
                               pattern = "^[^.].*\\.fa$", full.names = TRUE)
  if (length(alignment.files) == 0) stop("No expanded X5 alignments were found.")
  settings = list(
    run.TrimAl = run.TrimAl, trimAl.path = trimAl.path,
    min.external.percent = min.external.percent,
    min.column.gap.percent = min.column.gap.percent,
    min.coverage.percent = min.coverage.percent,
    min.coverage.bp = min.coverage.bp,
    min.alignment.length = min.alignment.length,
    min.taxa.alignment = min.taxa.alignment,
    max.alignment.gap.percent = max.alignment.gap.percent,
    manifest = .x5FileManifest(alignment.files)
  )
  .x5StageSettings(output.directory, "trim", settings, overwrite)
  .x5InvalidateDownstream(output.directory, "trim", overwrite)
  trimmed.directory = file.path(output.directory, "3_tree-alignments")
  excluded.directory = file.path(output.directory, "excluded", "trimmed")
  work.root = file.path(output.directory, "logs", "trimming")
  dir.create(file.path(output.directory, "tables"), recursive = TRUE,
             showWarnings = FALSE)
  .x5ClearDirectory(trimmed.directory, overwrite)
  .x5ClearDirectory(excluded.directory, overwrite)
  .x5ClearDirectory(work.root, overwrite)

  rows = list()
  decision.rows = list()
  for (index in seq_along(alignment.files)) {
    locus.id = sub("\\.fa$", "", basename(alignment.files[index]))
    collected.file = file.path(collected.directory, paste0(locus.id, ".rds"))
    if (!file.exists(collected.file)) stop("Missing collection record for ", locus.id)
    collected = readRDS(collected.file)
    result.file = file.path(trimmed.directory, paste0(locus.id, ".rds"))
    if (!overwrite && file.exists(result.file)) {
      result = readRDS(result.file)
      if (result$pass) {
        .x5ValidateAlignment(.x5ReadAlignment(result$alignment.file, "fasta"),
                             locus.id, result$sequence.ids)
      }
      rows[[index]] = result$summary
      decision.rows[[index]] = result$copy.decisions
      next
    }

    alignment = .x5ReadAlignment(alignment.files[index], "fasta")
    trim = .x5TrimCopies(
      alignment = alignment, copy.map = collected$Copy_map,
      run.TrimAl = run.TrimAl, trimAl.path = trimAl.path,
      work.directory = file.path(work.root, locus.id),
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.coverage.percent = min.coverage.percent,
      min.coverage.bp = min.coverage.bp,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      max.alignment.gap.percent = max.alignment.gap.percent,
      quiet = quiet
    )
    output.file = ""
    if (!is.null(trim$alignment)) {
      output.file = if (trim$pass) {
        file.path(trimmed.directory, paste0(locus.id, ".fa"))
      } else file.path(excluded.directory, paste0(locus.id, ".fa"))
      .x5WriteFasta(trim$alignment, output.file)
    }
    biological = collected$Copy_map$Role == "copy"
    initial.samples = length(unique(collected$Copy_map$Sample[biological]))
    initial.copies = sum(biological)
    survivors = if (is.null(trim$alignment)) character() else names(trim$alignment)
    copy.decisions = trim$removed
    remaining = setdiff(collected$Copy_map$Sequence_id, copy.decisions$Sequence_id)
    if (length(remaining) > 0) {
      copy.decisions = rbind(copy.decisions, data.frame(
        Sequence_id = remaining, Stage = "tree_preparation",
        Reason = if (trim$pass) "retained_for_tree" else "locus_failed_qc",
        Present_bp = NA_integer_, Coverage_percent = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    copy.decisions$Decision = ifelse(
      copy.decisions$Reason == "retained_for_tree", "retained", "excluded"
    )
    copy.decisions$Target = collected$Target
    copy.decisions = copy.decisions[, c(
      "Sequence_id", "Target", "Stage", "Decision", "Reason",
      "Present_bp", "Coverage_percent"
    )]
    summary = data.frame(
      Locus_id = locus.id, Target = collected$Target,
      Initial_samples = initial.samples, Initial_copies = initial.copies,
      Trimmed_samples = trim$sample.count, Trimmed_copies = trim$copy.count,
      Alignment_length = trim$length, Gap_percent = trim$gap.percent,
      Pass = trim$pass, Reason = trim$reason,
      Copy_evidence_removed = initial.copies > initial.samples &&
        trim$copy.count <= trim$sample.count,
      TrimAl_used = run.TrimAl, Alignment_file = output.file,
      stringsAsFactors = FALSE
    )
    result = list(
      stage = "trim", complete = TRUE, Locus_id = locus.id,
      Target = collected$Target, pass = trim$pass, reason = trim$reason,
      sequence.ids = survivors, alignment.file = output.file,
      copy.evidence.removed = summary$Copy_evidence_removed,
      copy.decisions = copy.decisions, summary = summary
    )
    saveRDS(result, result.file)
    rows[[index]] = summary
    decision.rows[[index]] = copy.decisions
  }
  summary = do.call(rbind, rows)
  copy.decisions = do.call(rbind, decision.rows)
  utils::write.table(summary,
                     file.path(output.directory, "tables", "trimming-summary.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(copy.decisions,
                     file.path(output.directory, "tables", "copy-decisions.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  message(sum(summary$Pass), " of ", nrow(summary),
          " targets passed copy-aware trimming.")
  invisible(summary)
}
