test_that("target names keep an underscore that belongs to the name", {
  locus.names <- c("gene_1", "gene_2")
  contig.names <- c("gene_1", "gene_2", "gene_2_1")

  expect_equal(.baseTargetName(contig.names, locus.names),
               c("gene_1", "gene_2", "gene_2"))
})

test_that("keepLonger takes the longer sequence and adds the new names", {
  current <- Biostrings::DNAStringSet(c(a = "AAAAAAAAAA", b = "CCCCCCCCCCCCCCC"))
  candidate <- Biostrings::DNAStringSet(c(a = "AAAAAAAAAAAAAAAAAAAA", c = "GGGGG"))

  result <- .keepLonger(current, candidate)

  expect_equal(sort(names(result)), c("a", "b", "c"))
  expect_equal(as.integer(Biostrings::width(result))[match("a", names(result))], 20L)
  expect_equal(as.integer(Biostrings::width(result))[match("b", names(result))], 15L)
})

test_that("a target with two contigs is left alone unless multi.copy is longest", {
  locus.names <- c("locus1", "locus2")
  old <- Biostrings::DNAStringSet(c(locus1 = "AAAAA",
                                    locus2 = "CCCCC",
                                    locus2_1 = "GGGGG"))
  new <- Biostrings::DNAStringSet(c(locus1 = "AAAAAAAAAA",
                                    locus2 = "TTTTTTTTTTTTTTT"))

  kept <- .mergeBinnedAssembly(old, new, locus.names, multi.copy = "keep")
  # locus1 is extended, both locus2 copies survive untouched
  expect_equal(as.integer(Biostrings::width(kept))[match("locus1", names(kept))], 10L)
  expect_equal(sum(.baseTargetName(names(kept), locus.names) == "locus2"), 2L)
  expect_false(any(Biostrings::width(kept) == 15L))

  replaced <- .mergeBinnedAssembly(old, new, locus.names, multi.copy = "longest")
  # the longest locus2 copy is beaten, the other copy stays
  expect_true(any(Biostrings::width(replaced) == 15L))
  expect_equal(sum(.baseTargetName(names(replaced), locus.names) == "locus2"), 2L)
})

test_that("a sample sequence is the bait only when it covers enough target", {
  reference <- Biostrings::DNAStringSet(c(
    locus1 = paste(rep("ACGT", 10), collapse = ""),
    locus2 = paste(rep("TGCA", 10), collapse = ""),
    locus3 = paste(rep("GGTT", 10), collapse = "")
  ))
  own <- Biostrings::DNAStringSet(c(locus1 = paste(rep("ACGT", 9), collapse = "")))
  draft <- Biostrings::DNAStringSet(c(locus2 = "TGCATGCA"))

  bait <- .buildBaitTable(reference.seqs = reference, own.contigs = own,
                          draft.contigs = draft, target.names = names(reference),
                          bait.source = "hybrid", min.bait.coverage = 0.5)

  # locus1 covers 90 percent, locus2 covers 20 percent, locus3 has no sequence
  expect_equal(bait$table$source, c("contig", "reference", "reference"))
  expect_equal(bait$table$locus, names(reference))
  expect_true(all(grepl("^bait[0-9]{6}$", names(bait$seqs))))
  expect_equal(length(bait$seqs), nrow(bait$table))

  # a longer draft contig is preferred to a shorter target contig
  long.draft <- Biostrings::DNAStringSet(c(locus1 = paste(rep("ACGT", 10), collapse = "")))
  both <- .buildBaitTable(reference.seqs = reference, own.contigs = own,
                          draft.contigs = long.draft, target.names = "locus1",
                          bait.source = "hybrid", min.bait.coverage = 0.5)
  expect_equal(both$table$source, "draft")

  # a rescued seed is used when nothing else covers the target
  seed <- Biostrings::DNAStringSet(c(locus3 = paste(rep("GGTT", 10), collapse = "")))
  rescued <- .buildBaitTable(reference.seqs = reference, own.contigs = own,
                             draft.contigs = draft, seed.contigs = seed,
                             target.names = names(reference),
                             bait.source = "hybrid", min.bait.coverage = 0.5)
  expect_equal(rescued$table$source, c("contig", "reference", "rescue"))

  # bait.source = reference ignores every sample sequence
  forced <- .buildBaitTable(reference.seqs = reference, own.contigs = own,
                            draft.contigs = draft, target.names = names(reference),
                            bait.source = "reference", min.bait.coverage = 0.5)
  expect_equal(unique(forced$table$source), "reference")
})

test_that("longestPerTarget keeps one contig per target", {
  locus.names <- c("locus1", "locus2")
  contigs <- Biostrings::DNAStringSet(c(locus1 = "AAAAA",
                                        locus1_1 = "AAAAAAAAAA",
                                        locus2 = "CCCCC",
                                        other = "GGGGG"))

  result <- .longestPerTarget(contigs, locus.names)

  expect_equal(sort(names(result)), c("locus1", "locus2"))
  expect_equal(as.integer(Biostrings::width(result))[match("locus1", names(result))], 10L)
})

test_that("a dinucleotide repeat is flagged and a varied sequence is not", {
  seqs <- Biostrings::DNAStringSet(c(
    repeated = paste(rep("AT", 60), collapse = ""),
    varied   = "ACGTTGCAAGGCTTACGATCCGTAAGCTTGCAATCGGATCCAAGTTGCA"
  ))

  expect_equal(unname(.lowComplexity(seqs)), c(TRUE, FALSE))
})
