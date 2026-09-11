make_x5_mafft = function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  executable = file.path(directory, "mafft")
  writeLines(c(
    "#!/bin/sh",
    "add=''",
    "last=''",
    "while [ \"$#\" -gt 0 ]; do",
    "  if [ \"$1\" = '--add' ]; then shift; add=\"$1\"; fi",
    "  last=\"$1\"",
    "  shift",
    "done",
    "cat \"$last\"",
    "if [ -n \"$add\" ]; then cat \"$add\"; fi"
  ), executable)
  Sys.chmod(executable, "0755")
  executable
}

test_that("copy-aware trimming counts each biological sample once", {
  alignment = Biostrings::DNAStringSet(c(
    S1 = "ACGTR---", S2 = "ACGT----", S3 = "ACGT----",
    S4 = "ACGT----", S5 = "----ACGT"
  ))
  map = data.frame(
    Sequence_id = names(alignment),
    Sample = c("A", "A", "B", "C", "D"),
    Role = "copy", stringsAsFactors = FALSE
  )
  result = .x5TrimCopies(
    alignment, map, run.TrimAl = FALSE,
    min.external.percent = 50, min.column.gap.percent = 100,
    min.coverage.percent = 0, min.coverage.bp = 0,
    min.alignment.length = 4, min.taxa.alignment = 4,
    max.alignment.gap.percent = 100
  )
  expect_true(result$pass)
  expect_equal(result$sample.count, 4)
  expect_equal(result$length, 5)
})

test_that("eligible copy splits are invariant to tree rooting", {
  tips = paste0("S", seq_len(8))
  map = data.frame(
    Sequence_id = tips,
    Sample = rep(LETTERS[1:4], 2),
    Role = "copy", Copy_group = "", stringsAsFactors = FALSE
  )
  tree = ape::read.tree(text = paste0(
    "((S1,S2,S3,S4)99:0.1,(S5,S6,S7,S8)99:0.1);"
  ))
  rerooted = ape::root(tree, outgroup = "S1", resolve.root = TRUE)
  first = .x5TreeSplits(tree, map, min.taxa.alignment = 4)
  second = .x5TreeSplits(rerooted, map, min.taxa.alignment = 4)
  expect_equal(sum(first$Eligible), 1)
  expect_equal(sum(second$Eligible), 1)
  expect_equal(first$Split_id[first$Eligible], second$Split_id[second$Eligible])
})

test_that("a deep single-copy split is not a copy concern", {
  tips = paste0("S", seq_len(8))
  map = data.frame(
    Sequence_id = tips, Sample = LETTERS[1:8], Role = "copy",
    Copy_group = "", stringsAsFactors = FALSE
  )
  tree = ape::read.tree(text =
    "((S1,S2,S3,S4)99:0.5,(S5,S6,S7,S8)99:0.5);")
  splits = .x5TreeSplits(
    tree, map, min.taxa.alignment = 4,
    min.split.branch.length = 0.1, min.split.branch.ratio = 1
  )
  expect_false(any(splits$Eligible))
  expect_false(any(splits$Deep_split_review))
})

test_that("a failed IQ-TREE command cannot reuse a stale tree", {
  root = tempfile("x5-failed-tree-")
  output = file.path(root, "output")
  trimmed.directory = file.path(output, "3_tree-alignments")
  collected.directory = file.path(output, "1_collected")
  tree.directory = file.path(output, "4_trees", "failed.target")
  tool.directory = file.path(root, "tools")
  for (directory in c(trimmed.directory, collected.directory,
                      tree.directory, tool.directory)) {
    dir.create(directory, recursive = TRUE)
  }
  executable = file.path(tool.directory, "iqtree2")
  writeLines(c("#!/bin/sh", "exit 2"), executable)
  Sys.chmod(executable, "0755")

  alignment = Biostrings::DNAStringSet(c(
    Species_alpha = "AAAA", Species_beta = "AAAC",
    Species_gamma = "AACC", Species_delta = "ACCC"
  ))
  alignment.file = file.path(trimmed.directory, "failed.target.fa")
  .x5WriteFasta(alignment, alignment.file)
  map = data.frame(
    Sequence_id = names(alignment), Target = "failed.target",
    Sample = LETTERS[1:4], Role = "copy", stringsAsFactors = FALSE
  )
  saveRDS(list(
    Locus_id = "failed.target", Target = "failed.target", Copy_map = map
  ), file.path(collected.directory, "failed.target.rds"))
  saveRDS(list(
    Locus_id = "failed.target", Target = "failed.target", pass = TRUE,
    alignment.file = alignment.file, sequence.ids = names(alignment)
  ), file.path(trimmed.directory, "failed.target.rds"))
  stale.tree = file.path(tree.directory, "failed.target.treefile")
  writeLines("((Species_alpha,Species_beta),(Species_gamma,Species_delta));",
             stale.tree)

  expect_error(inferParalogTrees(
    trimmed.directory = trimmed.directory,
    collected.directory = collected.directory,
    output.directory = output, tree.model = "JC",
    bootstrap.replicates = 1000, threads = 1, memory = 1,
    iqtree.path = tool.directory
  ), "IQ-TREE failed")
  expect_false(file.exists(stale.tree))
})

test_that("alignment expansion supports unchanged and missing-base targets", {
  root = tempfile("x5-alignment-paths-")
  output = file.path(root, "output")
  collected.directory = file.path(output, "1_collected")
  tool.directory = file.path(root, "tools")
  dir.create(collected.directory, recursive = TRUE)
  make_x5_mafft(tool.directory)
  base = Biostrings::DNAStringSet(c(Species_alpha = "ACGT",
                                    Species_beta = "ACGA"))
  additions = Biostrings::DNAStringSet(c(Species_gamma = "ACGT",
                                         Species_delta = "ACGA"))
  saveRDS(list(
    Locus_id = "base.target", Target = "base.target", Has_base = TRUE,
    Base_alignment = base, Additions = Biostrings::DNAStringSet()
  ), file.path(collected.directory, "base.target.rds"))
  saveRDS(list(
    Locus_id = "new.target", Target = "new.target", Has_base = FALSE,
    Base_alignment = Biostrings::DNAStringSet(), Additions = additions
  ), file.path(collected.directory, "new.target.rds"))
  guide.file = file.path(root, "targets.fa")
  .x5WriteFasta(Biostrings::DNAStringSet(c(new.target = "ACGT")), guide.file)

  summary = alignParalogCopies(
    collected.directory = collected.directory, target.file = guide.file,
    output.directory = output, threads = 1, mafft.path = tool.directory
  )
  expect_equal(summary$Expansion_path,
               c("unchanged_base", "new_alignment"))
  expect_equal(summary$Sequence_count, c(2, 2))
})

test_that("the best within-sample copy is retained without a supported split", {
  root = tempfile("x5-remove-paralog-sample-")
  output = file.path(root, "output")
  collected.directory = file.path(output, "1_collected")
  expanded.directory = file.path(output, "2_expanded")
  trimmed.directory = file.path(output, "3_tree-alignments")
  tree.directory = file.path(output, "4_trees", "target.one")
  for (directory in c(collected.directory, expanded.directory,
                      trimmed.directory, tree.directory)) {
    dir.create(directory, recursive = TRUE)
  }
  dir.create(file.path(output, "tables"), recursive = TRUE)
  utils::write.table(
    data.frame(Target = "target.one"),
    file.path(output, "tables", "target-map.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  alignment = Biostrings::DNAStringSet(c(
    A_1 = "ACGTACGT", A_2 = "TGCATGCA", B = "ACGTACGA",
    C = "ACGTACGC", D = "ACGTACGG", E = "ACGTACCC"
  ))
  map = data.frame(
    Sequence_id = names(alignment), Target = "target.one",
    Sample = c("A", "A", "B", "C", "D", "E"),
    Candidate_rank = c(1, 2, 1, 1, 1, 1),
    Source_contig = names(alignment), Role = "copy",
    stringsAsFactors = FALSE
  )
  expanded.file = file.path(expanded.directory, "target.one.fa")
  trimmed.file = file.path(trimmed.directory, "target.one.fa")
  .x5WriteFasta(alignment, expanded.file)
  .x5WriteFasta(alignment, trimmed.file)
  saveRDS(list(
    Locus_id = "target.one", Target = "target.one", Copy_map = map
  ), file.path(collected.directory, "target.one.rds"))
  copy.decisions = data.frame(
    Sequence_id = names(alignment), Target = "target.one",
    Stage = "tree_preparation", Decision = "retained",
    Reason = "retained_for_tree", Present_bp = NA_integer_,
    Coverage_percent = NA_real_, stringsAsFactors = FALSE
  )
  trim.summary = data.frame(
    Trimmed_samples = 5, Trimmed_copies = 6, Alignment_length = 8
  )
  saveRDS(list(
    Locus_id = "target.one", Target = "target.one", pass = TRUE,
    alignment.file = trimmed.file, sequence.ids = names(alignment),
    copy.evidence.removed = FALSE, copy.decisions = copy.decisions,
    summary = trim.summary
  ), file.path(trimmed.directory, "target.one.rds"))
  tree.file = file.path(tree.directory, "target.one.treefile")
  writeLines(paste0(
    "(((A_1:0.01,A_2:0.01)80:0.01,(B:0.01,C:0.01)80:0.01)",
    "80:0.01,(D:0.01,E:0.01)80:0.01);"
  ), tree.file)
  saveRDS(list(
    status = "complete", tree.file = tree.file,
    tip.ids = names(alignment)
  ), file.path(tree.directory, "result.rds"))

  decisions = separateParalogCopies(
    collected.directory = collected.directory,
    expanded.directory = expanded.directory,
    trimmed.directory = trimmed.directory,
    tree.directory = file.path(output, "4_trees"),
    output.directory = output, run.TrimAl = FALSE,
    min.external.percent = 0, min.column.gap.percent = 100,
    min.coverage.percent = 0, min.coverage.bp = 1,
    min.alignment.length = 1, min.taxa.alignment = 4,
    max.alignment.gap.percent = 100
  )
  expect_equal(decisions$Outcome,
               "retained_after_copy_selection")
  final = .x5ReadAlignment(file.path(
    output, "trimmed_all-markers", "target.one.phy"
  ), "phylip")
  expect_setequal(names(final), c("A", "B", "C", "D", "E"))
  copy.table = .x5ReadTable(file.path(output, "tables", "copy-decisions.tsv"))
  selected = copy.table[
    copy.table$Stage == "final" & copy.table$Sequence_id == "A_1",
  ]
  expect_equal(selected$Decision, "exported")
  redundant = copy.table[
    copy.table$Stage == "final" & copy.table$Sequence_id == "A_2",
  ]
  expect_equal(redundant$Decision, "excluded")
  expect_equal(redundant$Reason, "redundant_within_sample_copy")
})

test_that("workflow X5 exports both groups from one supported split", {
  root = tempfile("workflow-x5-")
  dir.create(root)
  alignment.directory = file.path(root, "base")
  paralog.directory = file.path(root, "paralogs")
  primary.directory = file.path(root, "primary")
  candidate.directory = file.path(root, "logs")
  output.directory = file.path(root, "output")
  tool.directory = file.path(root, "tools")
  for (directory in c(alignment.directory, paralog.directory,
                      primary.directory, candidate.directory)) {
    dir.create(directory)
  }
  make_x5_mafft(tool.directory)

  samples = c(
    "Celsiella_revocata", "Chimerella_corleone",
    "Cochranella_granulosa", "Vitreorana_antisthenesi"
  )
  base = Biostrings::DNAStringSet(setNames(rep("ACGTACGT", 4), samples))
  .x5WritePhylip(base, file.path(alignment.directory, "target.one.phy"))
  candidate.rows = list()
  for (index in seq_len(4)) {
    sample = samples[index]
    candidates = data.frame(
      Sample = sample, Target = "target.one",
      Contig = paste0("contig", seq_len(2)),
      Source_contig = paste0(sample, "_source", seq_len(2)),
      Candidate_rank = seq_len(2), Passes_absolute_filters = TRUE,
      Selected_for_primary = c(TRUE, FALSE), Saved_as_paralog = TRUE,
      Competitive_with_best = c(FALSE, TRUE), stringsAsFactors = FALSE
    )
    utils::write.csv(
      candidates,
      file.path(candidate.directory, paste0(sample, "_target-candidates.csv")),
      row.names = FALSE
    )
    sequences = Biostrings::DNAStringSet(c("ACGTACGT", "TGCATGCA"))
    names(sequences) = paste0(
      "target.one_|_", sample, "_|_copy", sprintf("%02d", seq_len(2))
    )
    .x5WriteFasta(sequences, file.path(paralog.directory, paste0(sample, ".fa")))
  }
  guide.file = file.path(root, "targets.fa")
  .x5WriteFasta(Biostrings::DNAStringSet(c(target.one = "ACGTACGT")),
                guide.file)

  collected = collectParalogCopies(
    alignment.directory = alignment.directory,
    paralog.directory = paralog.directory,
    primary.directory = primary.directory,
    candidate.directory = candidate.directory,
    output.directory = output.directory
  )
  expect_equal(collected$Initial_copies, 8)
  expect_equal(collected$Additions, 4)
  record = readRDS(file.path(output.directory, "1_collected", "target.one.rds"))
  expect_setequal(
    record$Copy_map$Sequence_id,
    as.vector(outer(samples, seq_len(2), paste, sep = "_"))
  )

  alignParalogCopies(
    collected.directory = file.path(output.directory, "1_collected"),
    target.file = guide.file, output.directory = output.directory,
    threads = 1, mafft.path = tool.directory
  )
  trimParalogCopies(
    expanded.directory = file.path(output.directory, "2_expanded"),
    collected.directory = file.path(output.directory, "1_collected"),
    output.directory = output.directory, run.TrimAl = FALSE,
    min.external.percent = 0, min.column.gap.percent = 100,
    min.coverage.percent = 0, min.coverage.bp = 1,
    min.alignment.length = 1, min.taxa.alignment = 4,
    max.alignment.gap.percent = 100
  )

  tree.directory = file.path(output.directory, "4_trees", "target.one")
  dir.create(tree.directory, recursive = TRUE)
  tree.file = file.path(tree.directory, "target.one.treefile")
  writeLines(paste0(
    "((Celsiella_revocata_1:0.01,Chimerella_corleone_1:0.01,",
    "Cochranella_granulosa_1:0.01,Vitreorana_antisthenesi_1:0.01)99:0.1,",
    "(Celsiella_revocata_2:0.01,Chimerella_corleone_2:0.01,",
    "Cochranella_granulosa_2:0.01,Vitreorana_antisthenesi_2:0.01)99:0.1);"
  ), tree.file)
  saveRDS(list(
    stage = "trees", complete = TRUE, status = "complete",
    tree.file = tree.file,
    tip.ids = c(paste0(samples, "_1"), paste0(samples, "_2")),
    summary = data.frame()
  ), file.path(tree.directory, "result.rds"))

  separation.arguments = list(
    collected.directory = file.path(output.directory, "1_collected"),
    expanded.directory = file.path(output.directory, "2_expanded"),
    trimmed.directory = file.path(output.directory, "3_tree-alignments"),
    tree.directory = file.path(output.directory, "4_trees"),
    output.directory = output.directory, run.TrimAl = FALSE,
    min.external.percent = 0, min.column.gap.percent = 100,
    min.coverage.percent = 0, min.coverage.bp = 1,
    min.alignment.length = 1, min.taxa.alignment = 4,
    max.alignment.gap.percent = 100
  )
  decisions = do.call(separateParalogCopies, separation.arguments)
  expect_equal(decisions$Outcome, "split")
  expect_equal(decisions$Output_marker_count, 2)
  expect_true(all(file.exists(file.path(
    output.directory, "trimmed_all-markers",
    c("target.one_copyA.phy", "target.one_copyB.phy")
  ))))
  copy.a = .x5ReadAlignment(file.path(
    output.directory, "trimmed_all-markers", "target.one_copyA.phy"
  ), "phylip")
  copy.b = .x5ReadAlignment(file.path(
    output.directory, "trimmed_all-markers", "target.one_copyB.phy"
  ), "phylip")
  expect_true(all(grepl("_[12]$", c(names(copy.a), names(copy.b)))))
  split.evidence = .x5ReadTable(file.path(
    output.directory, "tables", "split-evidence.tsv"
  ))
  expect_true(any(grepl("Celsiella_revocata_1", split.evidence$Side_a) |
                  grepl("Celsiella_revocata_1", split.evidence$Side_b)))
  resumed = do.call(separateParalogCopies, separation.arguments)
  expect_equal(resumed, decisions)
})
