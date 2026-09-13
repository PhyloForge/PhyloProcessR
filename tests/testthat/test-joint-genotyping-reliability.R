# Regression tests for workflow X1 joint genotyping. These cover the reference,
# cohort, and resume logic that runs before any external GATK command, so they
# need no bwa/samtools/gatk on the PATH.

# Writes a tiny gzipped single-sample GVCF header for the given loci and sample.
.write_test_gvcf = function(path, sample, loci) {
  con = gzfile(path, "wt")
  contig.lines = paste0("##contig=<ID=", loci, ",length=100>")
  writeLines(c("##fileformat=VCFv4.2", contig.lines,
               paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                       "INFO", "FORMAT", sample), collapse = "\t")), con)
  close(con)
}

# Builds a minimal dataset with a reference and one GVCF per sample.
.setup_cohort = function(samples, loci = c("locus1", "locus2")) {
  root = tempfile(); dir.create(root)
  reference.dir = file.path(root, "reference"); dir.create(reference.dir)
  reference = file.path(reference.dir, "reference.fa")
  writeLines(unlist(lapply(loci, function(l) c(paste0(">", l), paste(rep("A", 100), collapse = "")))),
             reference)
  file.create(paste0(reference, ".fai"))
  file.create(file.path(reference.dir, "reference.dict"))
  hap = file.path(root, "haplotype-caller"); dir.create(hap)
  for (s in samples) {
    d = file.path(hap, s); dir.create(d)
    g = file.path(d, "gatk4-haplotype-caller.g.vcf.gz")
    .write_test_gvcf(g, s, loci)
    file.create(paste0(g, ".tbi"))
  }
  list(root = root, reference = reference, hap = hap, loci = loci)
}


test_that("VCF header parsing reads the sample and contigs without loading records", {
  path = tempfile(fileext = ".vcf.gz")
  .write_test_gvcf(path, "SampleA", c("locus1", "locus10"))
  header = .vcfHeaderLines(path)
  expect_equal(.vcfSamples(header), "SampleA")
  expect_setequal(.vcfContigs(header), c("locus1", "locus10"))
})


test_that("a missing GVCF for an expected sample stops the run", {
  d = .setup_cohort(c("A", "B"))
  file.remove(file.path(d$hap, "B", "gatk4-haplotype-caller.g.vcf.gz"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap,
                    output.directory = file.path(d$root, "out"),
                    reference.path = d$reference,
                    sample.names = c("A", "B")),
    "Missing")
})


test_that("duplicate sample IDs across GVCFs stop the run", {
  d = .setup_cohort(c("A", "B"))
  # Rewrite B's GVCF to report the same internal sample ID as A
  .write_test_gvcf(file.path(d$hap, "B", "gatk4-haplotype-caller.g.vcf.gz"), "A", d$loci)
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap,
                    output.directory = file.path(d$root, "out"),
                    reference.path = d$reference,
                    sample.names = c("A", "B")),
    "Duplicate sample IDs")
})


test_that("a GVCF missing reference loci is rejected as a different reference", {
  d = .setup_cohort(c("A"))
  .write_test_gvcf(file.path(d$hap, "A", "gatk4-haplotype-caller.g.vcf.gz"), "A", c("other"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap,
                    output.directory = file.path(d$root, "out"),
                    reference.path = d$reference,
                    sample.names = "A"),
    "different reference")
})


test_that("disabling every final product is rejected before work begins", {
  d = .setup_cohort(c("A"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap,
                    output.directory = file.path(d$root, "out"),
                    reference.path = d$reference, sample.names = "A",
                    save.SNPs = FALSE, save.indels = FALSE, save.combined = FALSE),
    "save.SNPs")
})


test_that("a missing reference or its sidecars stops the run", {
  d = .setup_cohort(c("A"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap,
                    output.directory = file.path(d$root, "out"),
                    reference.path = file.path(d$root, "nope.fa"), sample.names = "A"),
    "Reference not found")
  file.remove(paste0(d$reference, ".fai"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap,
                    output.directory = file.path(d$root, "out"),
                    reference.path = d$reference, sample.names = "A"),
    "\\.fai")
})


test_that("a changed cohort refuses to resume without overwrite", {
  d = .setup_cohort(c("A", "B"))
  out = file.path(d$root, "out"); dir.create(out)
  # A record describing a different cohort than the current request
  saveRDS(list(samples = "different"), file.path(out, "cohort-record.rds"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap, output.directory = out,
                    reference.path = d$reference, sample.names = c("A", "B")),
    "changed")
})


test_that("filter thresholds default to usable numbers and reject non-finite values", {
  expect_equal(formals(jointGenotyping)$custom.SNP.QD, 2)
  expect_equal(eval(formals(jointGenotyping)$custom.SNP.MQRankSum), -12.5)
  d = .setup_cohort(c("A"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap, output.directory = file.path(d$root, "out"),
                    reference.path = d$reference, sample.names = "A", custom.SNP.QD = NA),
    "single finite number")
})


test_that("batch.size must be a positive integer", {
  d = .setup_cohort(c("A"))
  expect_error(
    jointGenotyping(haplotype.caller.directory = d$hap, output.directory = file.path(d$root, "out"),
                    reference.path = d$reference, sample.names = "A", batch.size = 0),
    "positive integer")
})


test_that("a fully completed cohort writes the reused sample-name map and skips work", {
  d = .setup_cohort(c("A", "B"))
  out = file.path(d$root, "out")
  # Pre-create every requested final product and its completion marker so no locus
  # is pending and no GATK command runs.
  for (sub in c("filtered-all", "filtered-snps", "filtered-indels")) {
    dir.create(file.path(out, sub), recursive = TRUE)
    for (locus in d$loci) writeLines("##fileformat=VCFv4.2", file.path(out, sub, paste0(locus, ".vcf")))
  }
  for (locus in d$loci) {
    .ensureDirectory(file.path(out, "completion", locus))
    .markStageComplete(file.path(out, "completion", locus), "jointGenotyping", "cohort=2")
  }
  expect_silent(
    jointGenotyping(haplotype.caller.directory = d$hap, output.directory = out,
                    reference.path = d$reference, sample.names = c("A", "B"), quiet = TRUE))
  map = file.path(out, "cohort-sample-map.txt")
  expect_true(file.exists(map))
  rows = read.table(map, sep = "\t", stringsAsFactors = FALSE)
  expect_equal(nrow(rows), 2)
  expect_setequal(rows[[1]], c("A", "B"))
})


# Writes a relaxed sequential phylip alignment.
.write_test_phylip = function(path, sequences) {
  writeLines(c(paste(length(sequences), nchar(sequences[[1]])),
               paste(names(sequences), unname(sequences))), path)
}


test_that("consensus mode needs supported phylip files and reads them all", {
  root = tempfile(); dir.create(root)
  writeLines("not an alignment", file.path(root, "report.txt"))
  expect_error(
    buildReference(reference.path = file.path(root, "reference", "reference.fa"),
                   reference.mode = "consensus", alignment.directory = root),
    "No phylip alignment files")

  .write_test_phylip(file.path(root, "geneA.phy"), c(t1 = "ACGT", t2 = "ACGA"))
  writeLines("garbage", file.path(root, "geneB.phy"))
  expect_error(
    buildReference(reference.path = file.path(root, "reference", "reference.fa"),
                   reference.mode = "consensus", alignment.directory = root),
    "Could not read alignment")
})


test_that("consensus mode rejects duplicate locus names and all-missing consensus", {
  # Same stem with two supported extensions collapses to one locus ID
  dup = tempfile(); dir.create(dup)
  .write_test_phylip(file.path(dup, "geneA.phy"), c(t1 = "ACGT", t2 = "ACGT"))
  .write_test_phylip(file.path(dup, "geneA.phylip"), c(t1 = "ACGT", t2 = "ACGT"))
  expect_error(
    buildReference(reference.path = file.path(dup, "reference", "reference.fa"),
                   reference.mode = "consensus", alignment.directory = dup),
    "Duplicate locus names")

  # An all-gap alignment produces an empty consensus
  empty = tempfile(); dir.create(empty)
  .write_test_phylip(file.path(empty, "geneA.phy"), c(t1 = "----", t2 = "----"))
  expect_error(
    buildReference(reference.path = file.path(empty, "reference", "reference.fa"),
                   reference.mode = "consensus", alignment.directory = empty),
    "empty or all-missing")
})


test_that("buildReference reuses an unchanged reference and refuses a changed one", {
  root = tempfile(); dir.create(root)
  source = file.path(root, "source.fa")
  writeLines(c(">locus1", paste(rep("A", 50), collapse = "")), source)
  reference.dir = file.path(root, "reference"); dir.create(reference.dir)
  reference = file.path(reference.dir, "reference.fa")
  outputs = c(reference, paste0(reference, c(".fai", ".amb", ".ann", ".bwt", ".pac", ".sa")),
              file.path(reference.dir, "reference.dict"))
  for (f in outputs) writeLines("x", f)

  # Record matching the current input reuses without touching any tool
  saveRDS(list(mode = "user", inputs = tools::md5sum(source)),
          file.path(reference.dir, "reference-build.rds"))
  expect_silent(buildReference(reference.path = reference, reference.mode = "user",
                               reference.file = source))

  # Changing the source content refuses reuse
  writeLines(c(">locus1", paste(rep("C", 50), collapse = "")), source)
  expect_error(buildReference(reference.path = reference, reference.mode = "user",
                              reference.file = source),
               "inputs changed")
})
