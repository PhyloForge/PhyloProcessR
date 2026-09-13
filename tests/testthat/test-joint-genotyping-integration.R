# End-to-end integration test for workflow X1. It runs the real GATK, BWA, and
# samtools pipeline on a tiny synthetic dataset, so it is skipped unless all three
# tools are on the PATH. Run it inside the PhyloProcessR conda environment.

random.seq = function(n) paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
revcomp = function(s) as.character(Biostrings::reverseComplement(Biostrings::DNAString(s)))

write.gz = function(lines, path) {
  con = gzfile(path, "wt"); writeLines(lines, con); close(con)
}

# Writes a relaxed sequential phylip alignment.
write.phylip = function(path, seqs) {
  writeLines(c(paste(length(seqs), nchar(seqs[[1]])),
               paste(names(seqs), unname(seqs))), path)
}

# Tiles paired reads across each sequence and writes one lane's R1/R2 files. The
# leading low-quality base forces FastqToSam to read Phred+33 rather than Solexa.
write.lane = function(dir, prefix, sequences, read.len = 100, frag = 250, step = 6) {
  high.chars = strsplit("FGHI", "")[[1]]
  make.qual = function() paste0("#", paste(sample(high.chars, read.len - 1, replace = TRUE), collapse = ""))
  r1 = character(0); r2 = character(0); n = 0
  for (locus in names(sequences)) {
    s = sequences[[locus]]
    for (i in seq(1, nchar(s) - frag, by = step)) {
      n = n + 1
      name = sprintf("M00001:1:FLOWCELL1:1:1101:%d:%d", i, n)
      read1 = substr(s, i, i + read.len - 1)
      read2 = revcomp(substr(s, i + frag - read.len, i + frag - 1))
      r1 = c(r1, paste0("@", name, " 1:N:0:1"), read1, "+", make.qual())
      r2 = c(r2, paste0("@", name, " 2:N:0:1"), read2, "+", make.qual())
    }
  }
  write.gz(r1, file.path(dir, paste0(prefix, "_R1.fastq.gz")))
  write.gz(r2, file.path(dir, paste0(prefix, "_R2.fastq.gz")))
}

vcf.records = function(path) {
  lines = readLines(path)
  lines[!startsWith(lines, "#")]
}
vcf.samples = function(path) {
  chrom = readLines(path); chrom = chrom[startsWith(chrom, "#CHROM")]
  tail(strsplit(chrom, "\t")[[1]], 2)
}


test_that("workflow X1 runs end to end against real GATK tools", {
  skip_if(Sys.which("gatk") == "", "gatk is not installed")
  skip_if(Sys.which("bwa") == "", "bwa is not installed")
  skip_if(Sys.which("samtools") == "", "samtools is not installed")

  set.seed(1)
  work = tempfile("x1-integration-")
  dir.create(work)
  old.wd = setwd(work)
  on.exit(setwd(old.wd), add = TRUE)

  # Reference loci and consensus alignment inputs (four identical copies).
  loci = c(locus1 = random.seq(600), locus2 = random.seq(600))
  dir.create("alignments")
  for (l in names(loci)) {
    write.phylip(file.path("alignments", paste0(l, ".phy")),
                 setNames(rep(loci[[l]], 4), paste0("t", 1:4)))
  }

  # sample1 carries a homozygous SNP at locus1 position 300; sample2 matches the
  # reference; locus2 is identical in both samples (a variant-free locus).
  snp.pos = 300
  alt = c(A = "T", C = "G", G = "C", T = "A")[substr(loci[["locus1"]], snp.pos, snp.pos)]
  s1.locus1 = loci[["locus1"]]; substr(s1.locus1, snp.pos, snp.pos) = alt
  sample1 = list(locus1 = s1.locus1, locus2 = loci[["locus2"]])
  sample2 = list(locus1 = loci[["locus1"]], locus2 = loci[["locus2"]])

  # sample1 is multilane; sample2 is single lane.
  dir.create("reads/sample1", recursive = TRUE)
  dir.create("reads/sample2", recursive = TRUE)
  write.lane("reads/sample1", "sample1_L001", sample1)
  write.lane("reads/sample1", "sample1_L002", sample1)
  write.lane("reads/sample2", "sample2_L001", sample2)

  mapping.dir = "sample-mapping"
  hap.dir = "haplotype-caller"
  geno.dir = "genotype-database"
  reference = file.path("reference", "reference.fa")

  prepareBAM(read.directory = "reads", output.directory = mapping.dir,
             threads = 2, memory = 4, quiet = TRUE)
  expected = list.dirs(mapping.dir, recursive = FALSE, full.names = FALSE)
  expect_setequal(expected, c("sample1", "sample2"))

  mapReferenceConsensus(mapping.directory = mapping.dir, alignment.directory = "alignments",
                        threads = 2, memory = 4, quiet = TRUE,
                        reference.path = reference, reference.mode = "consensus",
                        sample.names = expected)
  expect_length(Biostrings::readDNAStringSet(reference), 2)
  expect_true(all(file.exists(c(paste0(reference, ".fai"),
                                sub("\\.fa$", ".dict", reference)))))

  haplotypeCaller(mapping.directory = mapping.dir, output.directory = hap.dir,
                  reference.type = "consensus", reference.path = reference,
                  threads = 2, memory = 4, quiet = TRUE, sample.names = expected)

  jointGenotyping(haplotype.caller.directory = hap.dir, output.directory = geno.dir,
                  reference.path = reference, sample.names = expected,
                  threads = 2, memory = 4, quiet = TRUE)

  # The reused sample-name map lists the whole cohort.
  map.rows = read.table(file.path(geno.dir, "cohort-sample-map.txt"), sep = "\t",
                        stringsAsFactors = FALSE)
  expect_setequal(map.rows[[1]], c("sample1", "sample2"))

  # locus1 recovers the injected SNP with both sample columns; locus2 is empty.
  locus1.vcf = file.path(geno.dir, "filtered-snps", "locus1.vcf")
  locus2.vcf = file.path(geno.dir, "filtered-snps", "locus2.vcf")
  expect_true(file.exists(paste0(locus1.vcf, ".idx")))
  records1 = vcf.records(locus1.vcf)
  expect_gte(length(records1), 1)
  expect_equal(as.integer(strsplit(records1[1], "\t")[[1]][2]), snp.pos)
  expect_setequal(vcf.samples(locus1.vcf), c("sample1", "sample2"))
  expect_length(vcf.records(locus2.vcf), 0)

  # Resume with nothing changed skips completed loci.
  mtime.before = file.info(locus1.vcf)$mtime
  jointGenotyping(haplotype.caller.directory = hap.dir, output.directory = geno.dir,
                  reference.path = reference, sample.names = expected,
                  threads = 2, memory = 4, quiet = TRUE)
  expect_identical(file.info(locus1.vcf)$mtime, mtime.before)

  # An interrupted locus is rebuilt on the next run.
  unlink(file.path(geno.dir, "completion", "locus1"), recursive = TRUE)
  file.remove(locus1.vcf)
  jointGenotyping(haplotype.caller.directory = hap.dir, output.directory = geno.dir,
                  reference.path = reference, sample.names = expected,
                  threads = 2, memory = 4, quiet = TRUE)
  expect_true(file.exists(locus1.vcf))
})
