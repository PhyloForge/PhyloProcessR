#' Assess and separate putative paralog copies
#'
#' Evaluates supported unrooted tree splits, exports retained single-copy
#' markers, and exports both groups from one eligible two-copy split.
#'
#' @param collected.directory Directory made by [collectParalogCopies()].
#' @param expanded.directory Directory made by [alignParalogCopies()].
#' @param trimmed.directory Directory made by [trimParalogCopies()].
#' @param tree.directory Directory made by [inferParalogTrees()].
#' @param output.directory Workflow X5 output directory.
#' @param min.branch.support Inclusive ultrafast-bootstrap support threshold.
#' @param min.shared.samples Minimum samples represented on both split sides.
#' @param min.shared.sample.fraction Minimum shared fraction of the smaller side.
#' @param min.split.branch.length Branch-length threshold for review flags.
#' @param min.split.branch.ratio Branch-ratio threshold for review flags.
#' @param review.action Action for review-only single-copy targets.
#' @param run.TrimAl Run TrimAl during independent group trimming.
#' @param trimAl.path Directory that contains TrimAl, or NULL for PATH.
#' @param min.external.percent Minimum unique-sample occupancy at alignment edges.
#' @param min.column.gap.percent Missing-sample percentage that removes a column.
#' @param min.coverage.percent Minimum present-base percentage for each copy.
#' @param min.coverage.bp Minimum present-base count for each copy.
#' @param min.alignment.length Minimum retained alignment width.
#' @param min.taxa.alignment Minimum retained biological sample count.
#' @param max.alignment.gap.percent Maximum missing percentage in biological copies.
#' @param quiet Suppress external-tool console output when TRUE.
#' @param overwrite Replace X5 separation outputs when TRUE.
#'
#' @return A locus-decision summary, invisibly.
#' @export
separateParalogCopies = function(
    collected.directory = "data-analysis/paralog-analysis/1_collected",
    expanded.directory = "data-analysis/paralog-analysis/2_expanded",
    trimmed.directory = "data-analysis/paralog-analysis/3_tree-alignments",
    tree.directory = "data-analysis/paralog-analysis/4_trees",
    output.directory = "data-analysis/paralog-analysis",
    min.branch.support = 95,
    min.shared.samples = 3,
    min.shared.sample.fraction = 0.5,
    min.split.branch.length = 0.05,
    min.split.branch.ratio = 10,
    review.action = c("exclude", "retain"),
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

  review.action = match.arg(review.action)
  if (!is.finite(min.branch.support) || min.branch.support < 0 ||
      min.branch.support > 100) stop("min.branch.support must be from 0 to 100.")
  if (!is.finite(min.shared.sample.fraction) ||
      min.shared.sample.fraction < 0 || min.shared.sample.fraction > 1) {
    stop("min.shared.sample.fraction must be from 0 to 1.")
  }
  collected.files = list.files(collected.directory,
                                pattern = "^[^.].*\\.rds$", full.names = TRUE)
  if (length(collected.files) == 0) stop("No collected X5 targets were found.")
  dependency.files = c(
    collected.files,
    list.files(trimmed.directory, pattern = "^[^.].*\\.rds$", full.names = TRUE),
    list.files(tree.directory, pattern = "result\\.rds$", recursive = TRUE,
               full.names = TRUE),
    list.files(expanded.directory, pattern = "^[^.].*\\.fa$", full.names = TRUE),
    list.files(trimmed.directory, pattern = "^[^.].*\\.fa$", full.names = TRUE),
    list.files(tree.directory, pattern = "\\.treefile$", recursive = TRUE,
               full.names = TRUE)
  )
  settings = list(
    min.branch.support = min.branch.support,
    min.shared.samples = min.shared.samples,
    min.shared.sample.fraction = min.shared.sample.fraction,
    min.split.branch.length = min.split.branch.length,
    min.split.branch.ratio = min.split.branch.ratio,
    review.action = review.action, run.TrimAl = run.TrimAl,
    trimAl.path = trimAl.path,
    trim = list(
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.coverage.percent = min.coverage.percent,
      min.coverage.bp = min.coverage.bp,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      max.alignment.gap.percent = max.alignment.gap.percent
    ),
    manifest = .x5FileManifest(dependency.files)
  )
  .x5StageSettings(output.directory, "separate", settings, overwrite)
  .x5InvalidateDownstream(output.directory, "separate", overwrite)
  accepted.directory = file.path(output.directory, "trimmed_all-markers")
  fasta.directory = file.path(output.directory, "retained-fasta")
  excluded.untrimmed = file.path(output.directory, "excluded", "untrimmed")
  excluded.trimmed = file.path(output.directory, "excluded", "trimmed")
  provisional.directory = file.path(output.directory, "provisional-groups")
  result.directory = file.path(output.directory, "5_separated")
  work.root = file.path(output.directory, "logs", "group-trimming")
  dir.create(file.path(output.directory, "tables"), recursive = TRUE,
             showWarnings = FALSE)
  for (directory in c(accepted.directory, fasta.directory, excluded.untrimmed,
                      provisional.directory,
                      result.directory, work.root)) {
    .x5ClearDirectory(directory, overwrite)
  }
  dir.create(excluded.trimmed, recursive = TRUE, showWarnings = FALSE)

  target.map = .x5ReadTable(file.path(output.directory, "tables", "target-map.tsv"))
  all.products = c(target.map$Target, paste0(target.map$Target, "_copyA"),
                   paste0(target.map$Target, "_copyB"))
  collision.names = unique(all.products[duplicated(all.products) |
                                         duplicated(all.products, fromLast = TRUE)])
  safe.output.base = function(target, locus.id) {
    products = c(target, paste0(target, "_copyA"), paste0(target, "_copyB"))
    if (any(products %in% collision.names)) locus.id else .x5SafeTargetName(target)
  }

  trim.group = function(alignment, map, locus.id, label, use.trimal) {
    .x5TrimCopies(
      alignment = alignment, copy.map = map, run.TrimAl = use.trimal,
      trimAl.path = trimAl.path,
      work.directory = file.path(work.root, locus.id, label),
      min.external.percent = min.external.percent,
      min.column.gap.percent = min.column.gap.percent,
      min.coverage.percent = min.coverage.percent,
      min.coverage.bp = min.coverage.bp,
      min.alignment.length = min.alignment.length,
      min.taxa.alignment = min.taxa.alignment,
      max.alignment.gap.percent = max.alignment.gap.percent,
      quiet = quiet
    )
  }
  biological.names = function(alignment, map) {
    map = map[match(names(alignment), map$Sequence_id), , drop = FALSE]
    biological = map$Role == "copy"
    alignment = alignment[biological]
    map = map[biological, , drop = FALSE]
    if (anyDuplicated(map$Sample)) {
      stop("An exported marker would contain more than one copy for a sample.")
    }
    alignment
  }
  publish = function(alignment, marker) {
    fasta.file = file.path(fasta.directory, paste0(marker, ".fa"))
    phylip.file = file.path(accepted.directory, paste0(marker, ".phy"))
    .x5WriteFasta(alignment, fasta.file)
    .x5WritePhylip(alignment, phylip.file)
    list(fasta = normalizePath(fasta.file), phylip = normalizePath(phylip.file))
  }

  locus.rows = list()
  split.rows = list()
  membership.rows = list()
  output.rows = list()
  copy.rows = list()
  terminal.rows = list()

  for (index in seq_along(collected.files)) {
    collected = readRDS(collected.files[index])
    locus.id = unname(collected$Locus_id)
    result.file = file.path(result.directory, paste0(locus.id, ".rds"))
    if (!overwrite && file.exists(result.file)) {
      saved = readRDS(result.file)
      required.saved = c("locus.row", "split.rows", "membership.rows",
                         "output.rows", "copy.rows", "terminal.rows")
      complete.record = all(required.saved %in% names(saved))
      valid.outputs = complete.record &&
        (nrow(saved$output.rows) == 0 ||
         all(file.exists(saved$output.rows$Alignment_file)))
      valid.exclusion = complete.record &&
        (saved$locus.row$Outcome != "excluded" ||
         file.exists(saved$locus.row$Excluded_path))
      if (isTRUE(saved$complete) && complete.record &&
          valid.outputs && valid.exclusion) {
        locus.rows[[index]] = saved$locus.row
        if (nrow(saved$split.rows) > 0) {
          split.rows[[length(split.rows) + 1]] = saved$split.rows
        }
        if (nrow(saved$membership.rows) > 0) {
          membership.rows[[length(membership.rows) + 1]] = saved$membership.rows
        }
        if (nrow(saved$output.rows) > 0) {
          output.rows[[length(output.rows) + 1]] = saved$output.rows
        }
        copy.rows[[length(copy.rows) + 1]] = saved$copy.rows
        if (nrow(saved$terminal.rows) > 0) {
          terminal.rows[[length(terminal.rows) + 1]] = saved$terminal.rows
        }
        next
      }
    }
    split.start = length(split.rows)
    membership.start = length(membership.rows)
    output.start = length(output.rows)
    copy.start = length(copy.rows)
    terminal.start = length(terminal.rows)
    trim.file = file.path(trimmed.directory, paste0(locus.id, ".rds"))
    tree.result.file = file.path(tree.directory, locus.id, "result.rds")
    if (!file.exists(trim.file) || !file.exists(tree.result.file)) {
      stop("Completed trimming and tree records are required for ", locus.id)
    }
    trim = readRDS(trim.file)
    tree.result = readRDS(tree.result.file)
    expanded.file = file.path(expanded.directory, paste0(locus.id, ".fa"))
    if (!file.exists(expanded.file)) stop("Missing expanded alignment for ", locus.id)
    if (overwrite) {
      unlink(file.path(excluded.untrimmed, paste0(locus.id, ".fa")))
      if (trim$pass) {
        unlink(file.path(excluded.trimmed, paste0(locus.id, ".fa")))
      }
    }
    expanded = .x5ReadAlignment(expanded.file, "fasta")
    map = collected$Copy_map
    biological = map$Role == "copy"
    initial.samples = length(unique(map$Sample[biological]))
    initial.copies = sum(biological)
    outcome = "excluded"
    reason = "failed_qc"
    review.flag = FALSE
    eligible.count = 0L
    marker.count = 0L
    diagnostic.path = trim$alignment.file
    locus.splits = data.frame()
    exported.ids = character()
    redundant.copy.ids = character()

    if (trim$pass && identical(tree.result$status, "complete")) {
      tree = ape::read.tree(tree.result$tree.file)
      .x5ValidateTree(tree, trim$sequence.ids,
                      paste(collected$Target, "assessment tree"))
      terminal.edges = which(tree$edge[, 2] <= length(tree$tip.label))
      terminal.ids = tree$tip.label[tree$edge[terminal.edges, 2]]
      terminal.map = map[match(terminal.ids, map$Sequence_id), , drop = FALSE]
      terminal.rows[[length(terminal.rows) + 1]] = data.frame(
        Locus_id = locus.id, Target = collected$Target,
        Sequence_id = unname(terminal.ids),
        Sample = unname(terminal.map$Sample),
        Role = unname(terminal.map$Role),
        Branch_length = unname(tree$edge.length[terminal.edges]),
        stringsAsFactors = FALSE
      )
      survivor.map = map[map$Sequence_id %in% tree$tip.label, , drop = FALSE]
      locus.splits = .x5TreeSplits(
        tree, survivor.map, min.branch.support, min.shared.samples,
        min.shared.sample.fraction, min.split.branch.length,
        min.split.branch.ratio, min.taxa.alignment
      )
      if (nrow(locus.splits) > 0) {
        locus.splits$Locus_id = locus.id
        locus.splits$Target = collected$Target
        split.rows[[length(split.rows) + 1]] = locus.splits
        eligible.count = sum(locus.splits$Eligible)
        review.flag = any(locus.splits$Deep_split_review)
      }
      all.copy.samples = map$Sample[map$Role == "copy"]
      paralog.samples = unique(all.copy.samples[
        duplicated(all.copy.samples) |
          duplicated(all.copy.samples, fromLast = TRUE)
      ])
      reference.conflict = nrow(locus.splits) > 0 &&
        any(grepl("reference_conflict", locus.splits$Failure_reason,
                  fixed = TRUE))

      if (eligible.count == 1) {
        split = locus.splits[which(locus.splits$Eligible), , drop = FALSE]
        first.ids = split$Side_a[[1]]
        second.ids = split$Side_b[[1]]
        tuple.key = function(ids) {
          selected = survivor.map[survivor.map$Sequence_id %in% ids &
                                    survivor.map$Role == "copy", , drop = FALSE]
          paste(sort(paste(selected$Sample, selected$Candidate_rank,
                           selected$Source_contig, sep = "\r")), collapse = "\n")
        }
        if (tuple.key(first.ids) <= tuple.key(second.ids)) {
          groups = list(copyA = first.ids, copyB = second.ids)
        } else groups = list(copyA = second.ids, copyB = first.ids)
        group.results = list()
        for (label in names(groups)) {
          ids = groups[[label]]
          ids = ids[ids %in% trim$sequence.ids]
          group.map = map[map$Sequence_id %in% ids, , drop = FALSE]
          group.alignment = expanded[ids]
          group.results[[label]] = trim.group(
            group.alignment, group.map, locus.id, label, run.TrimAl
          )
          if (!is.null(group.results[[label]]$alignment)) {
            provisional = file.path(
              provisional.directory, paste0(locus.id, "_", label, ".fa")
            )
            .x5WriteFasta(group.results[[label]]$alignment, provisional)
          }
        }
        final.samples = lapply(group.results, function(result) {
          if (is.null(result$alignment)) return(character())
          group.map = map[match(names(result$alignment), map$Sequence_id), , drop = FALSE]
          unique(group.map$Sample[group.map$Role == "copy"])
        })
        shared = intersect(final.samples$copyA, final.samples$copyB)
        denominator = min(length(final.samples$copyA), length(final.samples$copyB))
        fraction = if (denominator == 0) 0 else length(shared) / denominator
        final.pass = all(vapply(group.results, function(result) result$pass,
                                logical(1))) &&
          length(shared) >= min.shared.samples &&
          fraction >= min.shared.sample.fraction
        for (label in names(groups)) {
          result = group.results[[label]]
          ids = if (is.null(result$alignment)) character() else names(result$alignment)
          membership.rows[[length(membership.rows) + 1]] = data.frame(
            Locus_id = locus.id, Target = collected$Target,
            Proposed_group = label, Sequence_id = groups[[label]],
            Survived_group_QC = groups[[label]] %in% ids,
            Exported = final.pass & groups[[label]] %in% ids,
            stringsAsFactors = FALSE
          )
        }
        if (final.pass) {
          base.name = safe.output.base(collected$Target, locus.id)
          for (label in names(groups)) {
            final = biological.names(group.results[[label]]$alignment, map)
            marker = paste0(base.name, "_", label)
            files = publish(final, marker)
            marker.count = marker.count + 1L
            ids = names(group.results[[label]]$alignment)
            exported.ids = c(exported.ids, ids)
            output.rows[[length(output.rows) + 1]] = data.frame(
              Output_marker = marker, Source_target = collected$Target,
              Copy_label = label, Alignment_file = files$phylip,
              Sample_count = length(final), Length = Biostrings::width(final)[1],
              Outcome = "split", Downstream_gene = "",
              Copy_identity_verified = TRUE, Unlinked_eligible = FALSE,
              stringsAsFactors = FALSE
            )
          }
          outcome = "split"
          reason = "supported_two_group_split"
        } else {
          reason = "split_failed_final_qc"
        }
      } else if (eligible.count > 1) {
        reason = "ambiguous_split"
      } else if (length(paralog.samples) > 0) {
        diagnostic = .x5ReadAlignment(trim$alignment.file, "fasta")
        diagnostic.map = map[
          match(names(diagnostic), map$Sequence_id), , drop = FALSE
        ]
        present.bp = rowSums(.x5PresentMatrix(diagnostic))
        selected.ids = character()
        for (sample in unique(diagnostic.map$Sample[
          diagnostic.map$Role == "copy"
        ])) {
          rows = which(diagnostic.map$Role == "copy" &
                         diagnostic.map$Sample == sample)
          ranks = diagnostic.map$Candidate_rank[rows]
          ranks[is.na(ranks)] = .Machine$integer.max
          order.rows = order(ranks, -present.bp[rows],
                             diagnostic.map$Source_contig[rows],
                             diagnostic.map$Sequence_id[rows])
          selected.ids = c(selected.ids,
                           diagnostic.map$Sequence_id[rows[order.rows[1]]])
        }
        redundant.copy.ids = diagnostic.map$Sequence_id[
          diagnostic.map$Role == "copy" &
            !(diagnostic.map$Sequence_id %in% selected.ids)
        ]
        clean.ids = c(
          diagnostic.map$Sequence_id[diagnostic.map$Role != "copy"],
          selected.ids
        )
        clean.map = map[map$Sequence_id %in% clean.ids, , drop = FALSE]
        clean = trim.group(
          diagnostic[clean.ids], clean.map, locus.id,
          "best-copies-selected", FALSE
        )
        if (clean$pass) {
          named = biological.names(clean$alignment, map)
          selected.map = map[match(names(named), map$Sequence_id), , drop = FALSE]
          names(named) = selected.map$Sample
          marker = safe.output.base(collected$Target, locus.id)
          files = publish(named, marker)
          marker.count = 1L
          exported.ids = names(clean$alignment)
          outcome = "retained_after_copy_selection"
          reason = "no_supported_deep_split_best_copy_selected"
          review.flag = FALSE
          diagnostic.path = files$fasta
          output.rows[[length(output.rows) + 1]] = data.frame(
            Output_marker = marker, Source_target = collected$Target,
            Copy_label = "", Alignment_file = files$phylip,
            Sample_count = length(named), Length = Biostrings::width(named)[1],
            Outcome = outcome, Downstream_gene = "",
            Copy_identity_verified = TRUE, Unlinked_eligible = TRUE,
            stringsAsFactors = FALSE
          )
        } else {
          reason = "failed_qc_after_copy_selection"
        }
      } else if (reference.conflict) {
        reason = "reference_conflict"
      } else if (review.flag && review.action == "exclude") {
        reason = "deep_split_review"
      } else {
        diagnostic = .x5ReadAlignment(trim$alignment.file, "fasta")
        final.map = map[map$Sequence_id %in% names(diagnostic), , drop = FALSE]
        final = trim.group(diagnostic, final.map, locus.id, "retained", FALSE)
        if (final$pass) {
          named = biological.names(final$alignment, map)
          marker = safe.output.base(collected$Target, locus.id)
          files = publish(named, marker)
          marker.count = 1L
          exported.ids = names(final$alignment)
          outcome = if (review.flag) "retained_with_review" else "retained"
          reason = if (review.flag) "deep_split_review" else "no_detected_copy_concern"
          output.rows[[length(output.rows) + 1]] = data.frame(
            Output_marker = marker, Source_target = collected$Target,
            Copy_label = "", Alignment_file = files$phylip,
            Sample_count = length(named), Length = Biostrings::width(named)[1],
            Outcome = outcome, Downstream_gene = "",
            Copy_identity_verified = TRUE, Unlinked_eligible = TRUE,
            stringsAsFactors = FALSE
          )
        } else reason = "failed_final_qc"
      }
    } else if (identical(tree.result$status, "insufficient_tree_information")) {
      reason = "insufficient_tree_information"
    } else if (identical(tree.result$status, "failed_tool")) {
      reason = "failed_tool"
    }

    if (outcome == "excluded") {
      file.copy(expanded.file,
                file.path(excluded.untrimmed, paste0(locus.id, ".fa")),
                overwrite = TRUE)
      excluded.trimmed.file = file.path(excluded.trimmed,
                                        paste0(locus.id, ".fa"))
      if (nzchar(diagnostic.path) && file.exists(diagnostic.path) &&
          normalizePath(diagnostic.path) !=
          normalizePath(excluded.trimmed.file, mustWork = FALSE)) {
        file.copy(diagnostic.path, excluded.trimmed.file, overwrite = TRUE)
      }
    }
    earlier.decisions = trim$copy.decisions
    earlier.decisions$Proposed_group = ""
    earlier.decisions$Exported_group = ""
    copy.rows[[length(copy.rows) + 1]] = earlier.decisions
    membership = if (length(membership.rows) == 0) data.frame() else
      do.call(rbind, membership.rows)
    for (sequence.id in map$Sequence_id) {
      proposed = ""
      exported.group = ""
      if (nrow(membership) > 0) {
        member = membership[membership$Locus_id == locus.id &
                              membership$Sequence_id == sequence.id, , drop = FALSE]
        if (nrow(member) == 1) {
          proposed = member$Proposed_group
          if (member$Exported) exported.group = proposed
        }
      }
      copy.rows[[length(copy.rows) + 1]] = data.frame(
        Sequence_id = sequence.id, Target = collected$Target,
        Stage = "final", Decision = if (sequence.id %in% exported.ids)
          "exported" else "excluded",
        Reason = if (sequence.id %in% exported.ids) outcome else
          if (sequence.id %in% redundant.copy.ids)
            "redundant_within_sample_copy" else reason,
        Present_bp = NA_integer_, Coverage_percent = NA_real_,
        Proposed_group = proposed, Exported_group = exported.group,
        stringsAsFactors = FALSE
      )
    }
    locus.rows[[index]] = data.frame(
      Locus_id = locus.id, Target = collected$Target,
      Initial_samples = initial.samples, Initial_copies = initial.copies,
      Trimmed_samples = trim$summary$Trimmed_samples,
      Trimmed_copies = trim$summary$Trimmed_copies,
      Trimmed_length = trim$summary$Alignment_length,
      Tree_status = tree.result$status, Deep_split_review = review.flag,
      Eligible_split_count = eligible.count, Outcome = outcome, Reason = reason,
      Output_marker_count = marker.count,
      Diagnostic_path = diagnostic.path,
      Excluded_path = if (outcome == "excluded")
        file.path(excluded.untrimmed, paste0(locus.id, ".fa")) else "",
      stringsAsFactors = FALSE
    )
    new.rows = function(rows, start) {
      if (length(rows) == start) return(data.frame())
      do.call(rbind, rows[seq.int(start + 1L, length(rows))])
    }
    saveRDS(list(
      stage = "separate", complete = reason != "failed_tool",
      locus.row = locus.rows[[index]],
      split.rows = new.rows(split.rows, split.start),
      membership.rows = new.rows(membership.rows, membership.start),
      output.rows = new.rows(output.rows, output.start),
      copy.rows = new.rows(copy.rows, copy.start),
      terminal.rows = new.rows(terminal.rows, terminal.start)
    ), result.file)
  }

  locus.decisions = do.call(rbind, locus.rows)
  split.evidence = if (length(split.rows) == 0) data.frame(
    Locus_id = character(), Target = character(), Split_id = character(),
    stringsAsFactors = FALSE) else do.call(rbind, split.rows)
  if (nrow(split.evidence) > 0) {
    split.evidence$Side_a = vapply(split.evidence$Side_a, paste,
                                   collapse = ";", character(1))
    split.evidence$Side_b = vapply(split.evidence$Side_b, paste,
                                   collapse = ";", character(1))
  }
  group.membership = if (length(membership.rows) == 0) data.frame(
    Locus_id = character(), Target = character(), Proposed_group = character(),
    Sequence_id = character(), Survived_group_QC = logical(),
    Exported = logical()) else do.call(rbind, membership.rows)
  output.markers = if (length(output.rows) == 0) data.frame(
    Output_marker = character(), Source_target = character(),
    Copy_label = character(), Alignment_file = character(),
    Sample_count = integer(), Length = integer(), Outcome = character(),
    Downstream_gene = character(), Copy_identity_verified = logical(),
    Unlinked_eligible = logical()) else do.call(rbind, output.rows)
  final.copy.decisions = do.call(rbind, copy.rows)
  terminal.branches = if (length(terminal.rows) == 0) data.frame(
    Locus_id = character(), Target = character(), Sequence_id = character(),
    Sample = character(), Role = character(), Branch_length = numeric()
  ) else do.call(rbind, terminal.rows)

  tables = file.path(output.directory, "tables")
  utils::write.table(locus.decisions, file.path(tables, "locus-decisions.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(split.evidence, file.path(tables, "split-evidence.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(group.membership, file.path(tables, "group-membership.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(output.markers, file.path(tables, "output-markers.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(final.copy.decisions,
                     file.path(tables, "copy-decisions.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(terminal.branches,
                     file.path(tables, "terminal-branches.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  message(sum(locus.decisions$Outcome %in% c(
    "retained", "retained_with_review",
    "retained_after_copy_selection", "split"
  )),
          " source targets produced ", nrow(output.markers),
          " accepted marker alignments. Tables: ", tables)
  invisible(locus.decisions)
}
