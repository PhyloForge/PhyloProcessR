#' Infer gene trees for workflow X5
#'
#' Runs IQ-TREE for each copy-aware alignment that has sufficient sample and
#' variable-site information. Command failures stop the stage and cannot reuse
#' an older tree.
#'
#' @param trimmed.directory Directory made by [trimParalogCopies()].
#' @param collected.directory Directory made by [collectParalogCopies()].
#' @param output.directory Workflow X5 output directory.
#' @param tree.model IQ-TREE model setting.
#' @param bootstrap.replicates Number of ultrafast bootstrap replicates.
#' @param tree.seed IQ-TREE random seed.
#' @param threads Number of IQ-TREE threads.
#' @param memory Memory limit in GB.
#' @param iqtree.path Directory that contains IQ-TREE, or NULL for PATH.
#' @param iqtree.executable Compatible IQ-TREE executable name.
#' @param quiet Suppress external-tool console output when TRUE.
#' @param overwrite Replace X5 tree outputs when TRUE.
#'
#' @return A tree-inference summary, invisibly.
#' @export
inferParalogTrees = function(
    trimmed.directory = "data-analysis/paralog-analysis/3_tree-alignments",
    collected.directory = "data-analysis/paralog-analysis/1_collected",
    output.directory = "data-analysis/paralog-analysis",
    tree.model = "MFP",
    bootstrap.replicates = 1000,
    tree.seed = 12345,
    threads = 8,
    memory = 40,
    iqtree.path = NULL,
    iqtree.executable = "iqtree2",
    quiet = TRUE,
    overwrite = FALSE) {

  .validateResources(threads = threads, memory = memory)
  if (bootstrap.replicates < 1 || bootstrap.replicates != as.integer(bootstrap.replicates)) {
    stop("bootstrap.replicates must be a positive integer.")
  }
  trim.records = list.files(trimmed.directory, pattern = "^[^.].*\\.rds$",
                            full.names = TRUE)
  if (length(trim.records) == 0) stop("No X5 trimming records were found.")
  trimmed.alignments = vapply(trim.records, function(file) {
    result = readRDS(file)
    if (is.null(result$alignment.file)) "" else result$alignment.file
  }, character(1))
  settings = list(
    tree.model = tree.model,
    bootstrap.replicates = as.integer(bootstrap.replicates),
    tree.seed = as.integer(tree.seed), threads = as.integer(threads),
    memory = memory, iqtree.path = iqtree.path,
    iqtree.executable = iqtree.executable,
    manifest = .x5FileManifest(c(trim.records, trimmed.alignments))
  )
  .x5StageSettings(output.directory, "trees", settings, overwrite)
  .x5InvalidateDownstream(output.directory, "trees", overwrite)
  tree.directory = file.path(output.directory, "4_trees")
  .x5ClearDirectory(tree.directory, overwrite)
  executable = NULL
  tool.version = NA_character_
  rows = list()
  write.run.info = function(status) {
    package.version = tryCatch(
      as.character(utils::packageVersion("PhyloProcessR")),
      error = function(error) "development source"
    )
    writeLines(c(
      paste0("status=", status),
      paste0("PhyloProcessR_version=", package.version),
      paste0("R_version=", R.version.string),
      paste0("IQ_TREE_version=", tool.version),
      paste0("tree_model=", tree.model),
      paste0("bootstrap_replicates=", as.integer(bootstrap.replicates)),
      paste0("tree_seed=", as.integer(tree.seed)),
      paste0("threads=", as.integer(threads)),
      paste0("memory_GB=", memory)
    ), file.path(output.directory, "run-info.txt"))
  }
  write.run.info("incomplete")

  for (index in seq_along(trim.records)) {
    trim = readRDS(trim.records[index])
    locus.directory = file.path(tree.directory, trim$Locus_id)
    dir.create(locus.directory, recursive = TRUE, showWarnings = FALSE)
    result.file = file.path(locus.directory, "result.rds")
    if (!overwrite && file.exists(result.file)) {
      result = readRDS(result.file)
      if (identical(result$status, "complete")) {
        tree = ape::read.tree(result$tree.file)
        .x5ValidateTree(tree, result$tip.ids,
                        paste(trim$Locus_id, "saved tree"))
      }
      rows[[index]] = result$summary
      next
    }

    if (!trim$pass) {
      summary = data.frame(
        Locus_id = trim$Locus_id, Target = trim$Target,
        Tree_status = "not_run_failed_qc", Tree_file = "",
        Tip_count = 0, stringsAsFactors = FALSE
      )
      saveRDS(list(stage = "trees", complete = TRUE,
                   status = "not_run_failed_qc", tree.file = "",
                   tip.ids = character(), summary = summary), result.file)
      rows[[index]] = summary
      next
    }
    collected = readRDS(file.path(collected.directory,
                                  paste0(trim$Locus_id, ".rds")))
    alignment = .x5ReadAlignment(trim$alignment.file, "fasta")
    if (!.x5VariableSite(alignment, collected$Copy_map)) {
      summary = data.frame(
        Locus_id = trim$Locus_id, Target = trim$Target,
        Tree_status = "insufficient_tree_information", Tree_file = "",
        Tip_count = length(alignment), stringsAsFactors = FALSE
      )
      saveRDS(list(stage = "trees", complete = TRUE,
                   status = "insufficient_tree_information", tree.file = "",
                   tip.ids = names(alignment), summary = summary), result.file)
      rows[[index]] = summary
      next
    }

    if (is.null(executable)) {
      executable = .toolCommand(iqtree.executable, iqtree.path)
      version.output = try(.runCommandOutput(
        paste(executable, "--version"), task = "IQ-TREE version check"
      ), silent = TRUE)
      if (!inherits(version.output, "try-error")) {
        tool.version = paste(version.output, collapse = " ")
      }
    }
    prefix = file.path(locus.directory, trim$Locus_id)
    stale = list.files(locus.directory,
                       pattern = paste0("^", trim$Locus_id, "\\."),
                       full.names = TRUE)
    if (length(stale) > 0) unlink(stale)
    stdout.log = file.path(locus.directory, "iqtree.stdout.log")
    stderr.log = file.path(locus.directory, "iqtree.stderr.log")
    command = paste(
      executable, "-s", shQuote(trim$alignment.file),
      "-m", shQuote(tree.model), "-B", as.integer(bootstrap.replicates),
      "-bnni -seed", as.integer(tree.seed), "-nt", as.integer(threads),
      "-mem", paste0(memory, "G"), "-keep-ident -pre", shQuote(prefix)
    )
    failed = try(.runCommand(
      paste(command, ">", shQuote(stdout.log)), quiet = quiet,
      task = paste(trim$Target, "IQ-TREE"), keep.stdout = TRUE,
      stderr.log = stderr.log
    ), silent = TRUE)
    tree.file = paste0(prefix, ".treefile")
    if (inherits(failed, "try-error")) {
      unlink(tree.file)
      summary = data.frame(
        Locus_id = trim$Locus_id, Target = trim$Target,
        Tree_status = "failed_tool", Tree_file = "",
        Tip_count = length(alignment), stringsAsFactors = FALSE
      )
      saveRDS(list(stage = "trees", complete = FALSE, status = "failed_tool",
                   tree.file = "", tip.ids = names(alignment),
                   error = as.character(failed), log.file = stderr.log,
                   summary = summary), result.file)
      write.run.info("failed_tool")
      stop("IQ-TREE failed for ", trim$Target, ". See ", stderr.log)
    }
    if (!file.exists(tree.file) || file.info(tree.file)$size == 0) {
      stop("IQ-TREE did not create a tree for ", trim$Target)
    }
    tree = try(ape::read.tree(tree.file), silent = TRUE)
    if (inherits(tree, "try-error")) {
      stop("IQ-TREE output has an invalid tip set for ", trim$Target)
    }
    .x5ValidateTree(tree, names(alignment), paste(trim$Target, "IQ-TREE output"))
    summary = data.frame(
      Locus_id = trim$Locus_id, Target = trim$Target,
      Tree_status = "complete", Tree_file = normalizePath(tree.file),
      Tip_count = length(tree$tip.label), stringsAsFactors = FALSE
    )
    result = list(
      stage = "trees", complete = TRUE, status = "complete",
      tree.file = normalizePath(tree.file), tip.ids = names(alignment),
      command = command, tool.version = tool.version, summary = summary
    )
    saveRDS(result, result.file)
    rows[[index]] = summary
  }
  summary = do.call(rbind, rows)
  utils::write.table(summary,
                     file.path(output.directory, "tables", "tree-summary.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  write.run.info("complete")
  message(sum(summary$Tree_status == "complete"), " trees inferred. Output: ",
          tree.directory)
  invisible(summary)
}
