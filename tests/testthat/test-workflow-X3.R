# Regression tests for the workflow X3 legacy-integration fixes. Tests assert the
# fixed behaviour (not the pre-fix faults recorded in the review plan). The
# addLegacyAlignments tests replace BLAST and MAFFT with small shell stubs, so no
# real tool or biological result is exercised.

# Write a PHYLIP alignment from a named character vector of equal-length strings.
write_test_phylip <- function(sequences, file) {
  alignment <- as.matrix(ape::as.DNAbin(
    strsplit(as.character(Biostrings::DNAStringSet(sequences)), "")))
  PhyloProcessR::writePhylip(alignment = alignment, file = file,
                             interleave = FALSE, strict = FALSE)
}

# blastn stub that writes a single outfmt-6 hit whose subject is a fixed target.
# makeblastdb is a no-op that exits 0. Both are named to match the command the
# function builds (blast.path + executable).
make_blast_stubs <- function(directory, subject) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("#!/bin/sh", "exit 0"), file.path(directory, "makeblastdb"))
  writeLines(c(
    "#!/bin/sh",
    "out=''; prev=''",
    "for arg in \"$@\"; do",
    "  if [ \"$prev\" = \"-out\" ]; then out=\"$arg\"; fi",
    "  prev=\"$arg\"",
    "done",
    paste0("printf 'q\\t", subject,
           "\\t100\\t100\\t0\\t0\\t1\\t100\\t1\\t100\\t0.0\\t200\\t100\\t100\\t0\\n' > \"$out\"")
  ), file.path(directory, "blastn"))
  Sys.chmod(file.path(directory, "makeblastdb"), "0755")
  Sys.chmod(file.path(directory, "blastn"), "0755")
  directory
}

# blastn stub that reports no hit (writes an empty output file).
make_blast_stubs_nohit <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  writeLines(c("#!/bin/sh", "exit 0"), file.path(directory, "makeblastdb"))
  writeLines(c(
    "#!/bin/sh",
    "out=''; prev=''",
    "for arg in \"$@\"; do",
    "  if [ \"$prev\" = \"-out\" ]; then out=\"$arg\"; fi",
    "  prev=\"$arg\"",
    "done",
    ": > \"$out\""
  ), file.path(directory, "blastn"))
  Sys.chmod(file.path(directory, "makeblastdb"), "0755")
  Sys.chmod(file.path(directory, "blastn"), "0755")
  directory
}

# mafft --add stub: concatenate the two input FASTA files (equal-width fixtures
# give an equal-width "alignment"), preserving the tagged names the function set.
make_mafft_stub <- function(directory) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "#!/bin/sh",
    "for arg in \"$@\"; do",
    "  case \"$arg\" in",
    "    *.fa) if [ -f \"$arg\" ]; then cat \"$arg\"; fi ;;",
    "  esac",
    "done"
  ), file.path(directory, "mafft"))
  Sys.chmod(file.path(directory, "mafft"), "0755")
  directory
}

write_nexus <- function(file, sets) {
  writeLines(c(
    "#NEXUS",
    "BEGIN DATA;",
    "  DIMENSIONS NTAX=5 NCHAR=12;",
    "  FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=. ;",
    "  MATRIX",
    "  t1 ACGTACGTACGT",
    "  t2 ...A...A...A",
    "  t3 ACGTACGTACGT",
    "  t4 ACGTACGTACGT",
    "  t5 ACGTACGTACGT",
    "  ;",
    "END;",
    "BEGIN SETS;",
    sets,
    "END;"
  ), file)
}


test_that(".keepMostInformativeRows keeps the best copy without emptying the set", {
  # Two same-named rows, the second more informative: the pre-fix code produced
  # align[-integer(0)] and emptied the alignment.
  pair <- Biostrings::DNAStringSet(c(A = "AC--", A = "ACGT", B = "AAAA"))
  kept <- PhyloProcessR:::.keepMostInformativeRows(pair)
  expect_equal(length(kept), 2L)
  expect_equal(as.character(kept)[names(kept) == "A"], c(A = "ACGT"))
  expect_true("B" %in% names(kept))

  # Three copies, first most informative, and ties keep the first occurrence.
  triple <- Biostrings::DNAStringSet(c(X = "ACGT", X = "AC--", X = "A---", Y = "AAAA"))
  kept3 <- PhyloProcessR:::.keepMostInformativeRows(triple)
  expect_equal(length(kept3), 2L)
  expect_equal(as.character(kept3)[names(kept3) == "X"], c(X = "ACGT"))
})


test_that("convertNexusPartitions parses steps, singletons, and rejects bad charsets", {
  root <- tempfile("x3-nexus-")
  dir.create(root)
  nexus <- file.path(root, "matrix.nex")
  write_nexus(nexus, c(
    "  charset codon1 = 1-12\\3 ;",
    "  charset g2 = 1-4 ;",
    "  charset g3 = 7 ;",
    "  charset bad = 99-200 ;"))
  out <- file.path(root, "loci")

  convertNexusPartitions(nexus.file = nexus, output.directory = out,
                         output.format = "phylip", min.taxa.alignment = 4,
                         max.missing.percent = 100, overwrite = TRUE, quiet = TRUE)

  # The out-of-range charset is rejected; the three valid ones are written.
  expect_setequal(sub("\\.phy$", "", list.files(out)), c("codon1", "g2", "g3"))

  codon1 <- as.character(Biostrings::DNAStringSet(
    Biostrings::readDNAMultipleAlignment(file.path(out, "codon1.phy"), format = "phylip")))
  # Columns 1,4,7,10 of ACGTACGTACGT are A,T,G,C.
  expect_equal(unname(codon1["t1"]), "ATGC")
  # MATCHCHAR '.' for t2 is expanded from t1 before extraction.
  expect_equal(unname(codon1["t2"]), "AAGC")

  g3 <- as.character(Biostrings::DNAStringSet(
    Biostrings::readDNAMultipleAlignment(file.path(out, "g3.phy"), format = "phylip")))
  expect_equal(unname(nchar(g3["t1"])), 1L)
  expect_equal(unname(g3["t1"]), "G")
})


test_that("convertNexusPartitions stops on duplicate charset names", {
  root <- tempfile("x3-nexus-dup-")
  dir.create(root)
  nexus <- file.path(root, "matrix.nex")
  write_nexus(nexus, c("  charset dup = 1-4 ;", "  charset dup = 5-8 ;"))
  expect_error(
    convertNexusPartitions(nexus.file = nexus,
                           output.directory = file.path(root, "loci"),
                           min.taxa.alignment = 4, overwrite = TRUE, quiet = TRUE),
    "Duplicate charset")
})


test_that("convertNexusPartitions keeps an alignment with exactly the minimum taxa", {
  root <- tempfile("x3-nexus-min-")
  dir.create(root)
  nexus <- file.path(root, "matrix.nex")
  # Drop one taxon with max.missing.percent so exactly four remain for g2.
  write_nexus(nexus, "  charset g2 = 1-4 ;")
  out <- file.path(root, "loci")
  # t2 is all matchchar-derived but still has bases after expansion; keep all 5.
  convertNexusPartitions(nexus.file = nexus, output.directory = out,
                         min.taxa.alignment = 5, overwrite = TRUE, quiet = TRUE)
  expect_true(file.exists(file.path(out, "g2.phy")))
})


test_that("gatherUnlinked preserves loci absent from the metadata", {
  root <- tempfile("x3-gather-")
  dir.create(root)
  genes <- file.path(root, "genes")
  exons <- file.path(root, "exons")
  out <- file.path(root, "unlinked")
  dir.create(genes); dir.create(exons)

  # Metadata maps only marker1/marker2 to geneA. locusX is absent from it.
  metadata <- file.path(root, "metadata.txt")
  write.table(data.frame(marker = c("marker1", "marker2"), gene = c("geneA", "geneA")),
              metadata, sep = "\t", row.names = FALSE, quote = FALSE)

  write_test_phylip(c(t1 = "ACGT", t2 = "ACGT"), file.path(genes, "geneA.phy"))
  write_test_phylip(c(t1 = "ACGT", t2 = "ACGT"), file.path(exons, "marker1.phy"))
  write_test_phylip(c(t1 = "ACGT", t2 = "ACGT"), file.path(exons, "marker2.phy"))
  write_test_phylip(c(t1 = "ACGT", t2 = "ACGT"), file.path(exons, "locusX.phy"))

  gatherUnlinked(gene.alignment.directory = genes, exon.alignment.directory = exons,
                 output.directory = out, feature.gene.names = metadata, overwrite = TRUE)

  produced <- sub("\\.phy$", "", list.files(out))
  # The concatenated gene and the unmatched locus are present; the exons folded
  # into geneA are not duplicated.
  expect_true("geneA" %in% produced)
  expect_true("locusX" %in% produced)
  expect_false("marker1" %in% produced)
})


test_that("gatherUnlinked stops when a gene and exon name collide", {
  root <- tempfile("x3-gather-clash-")
  dir.create(root)
  genes <- file.path(root, "genes")
  exons <- file.path(root, "exons")
  dir.create(genes); dir.create(exons)
  metadata <- file.path(root, "metadata.txt")
  write.table(data.frame(marker = "other", gene = "shared"),
              metadata, sep = "\t", row.names = FALSE, quote = FALSE)
  write_test_phylip(c(t1 = "ACGT", t2 = "ACGT"), file.path(genes, "shared.phy"))
  # An unmatched exon named "shared" would overwrite the gene output.
  write_test_phylip(c(t1 = "ACGT", t2 = "ACGT"), file.path(exons, "shared.phy"))
  expect_error(
    gatherUnlinked(gene.alignment.directory = genes, exon.alignment.directory = exons,
                   output.directory = file.path(root, "out"),
                   feature.gene.names = metadata, overwrite = TRUE),
    "collide")
})


test_that("addLegacyAlignments resolves the target by exact locus ID", {
  skip_on_os("windows")
  root <- tempfile("x3-exact-")
  dir.create(root)
  captures <- file.path(root, "captures")
  legacy <- file.path(root, "legacy")
  dir.create(captures); dir.create(legacy)

  # Two capture loci whose names share a prefix: a substring match would pick both.
  write_test_phylip(c(capA = "ACGTACGT", capB = "ACGTACGT"), file.path(captures, "locus1.phy"))
  write_test_phylip(c(capA = "TTTTTTTT", capB = "TTTTTTTT"), file.path(captures, "locus10.phy"))
  write_test_phylip(c(legZ = "ACGTACGT"), file.path(legacy, "legloc.phy"))

  target.fa <- file.path(root, "targets.fa")
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(c(locus1 = "ACGTACGT", locus10 = "TTTTTTTT")), target.fa)

  blast.tools <- make_blast_stubs(file.path(root, "blast"), subject = "locus1")
  mafft.tools <- make_mafft_stub(file.path(root, "mafft"))
  out.base <- file.path(root, "untrimmed_legacy")

  addLegacyAlignments(alignment.directory = captures, alignment.format = "phylip",
                      output.directory = out.base, legacy.directory = legacy,
                      legacy.format = "phylip", target.markers = target.fa,
                      merge = "None", overwrite = TRUE, quiet = TRUE,
                      blast.path = blast.tools, mafft.path = mafft.tools)

  # Exactly locus1 was resolved, not locus10.
  expect_true(file.exists(file.path(paste0(out.base, "-only"), "locus1.phy")))
  expect_false(file.exists(file.path(paste0(out.base, "-only"), "locus10.phy")))
  # The completion summary marks a finished run.
  expect_true(file.exists(paste0(out.base, "-integration_summary.txt")))
})


test_that("addLegacyAlignments keeps the best uncaptured duplicate and cleans scratch", {
  skip_on_os("windows")
  root <- tempfile("x3-uncap-")
  dir.create(root)
  captures <- file.path(root, "captures")
  legacy <- file.path(root, "legacy")
  dir.create(captures); dir.create(legacy)
  write_test_phylip(c(capA = "ACGTACGT", capB = "ACGTACGT"), file.path(captures, "locus1.phy"))
  # Legacy locus has a duplicate taxon; the second copy is more informative.
  write_test_phylip(c(dupT = "AC------", dupT = "ACGTACGT", other = "ACGTACGT"),
                    file.path(legacy, "orphan.phy"))

  target.fa <- file.path(root, "targets.fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(locus1 = "ACGTACGT")), target.fa)
  blast.tools <- make_blast_stubs_nohit(file.path(root, "blast"))
  out.base <- file.path(root, "untrimmed_legacy")

  owd <- getwd(); setwd(root); on.exit(setwd(owd), add = TRUE)
  addLegacyAlignments(alignment.directory = captures, alignment.format = "phylip",
                      output.directory = out.base, legacy.directory = legacy,
                      legacy.format = "phylip", target.markers = target.fa,
                      merge = "None", include.uncaptured.legacy = TRUE,
                      overwrite = TRUE, quiet = TRUE,
                      blast.path = blast.tools, mafft.path = NULL)

  orphan <- Biostrings::DNAStringSet(
    Biostrings::readDNAMultipleAlignment(file.path(paste0(out.base, "-only"), "orphan.phy"),
                                         format = "phylip"))
  # The duplicate is reduced to the informative copy, not emptied.
  expect_equal(sum(names(orphan) == "dupT"), 1L)
  expect_equal(unname(as.character(orphan)["dupT"]), "ACGTACGT")
  expect_true("other" %in% names(orphan))
  # No scratch files were left in the working directory.
  expect_length(list.files(root, pattern = "_query\\.fa$|_blast-match\\.txt$|blast_db"), 0L)
})
