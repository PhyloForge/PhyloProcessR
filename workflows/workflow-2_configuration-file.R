#################################################
## Configuration file for PhyloProcessR
#################################################

# Package version
#########################
# TRUE installs the latest development/beta version from GitHub before running.
# Keep FALSE for reproducible analyses that use the already installed version.
install.latest.github = FALSE

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
working.directory = "/PATH/TO/where/the/stuff/will/happen"
#The processed reads directory
processed.reads = "processed-reads"
#The read folder within processed reads to assemble, pe-merged-reads recommended
assembly.reads = "pe-merged-reads"
# The target markers reference file
target.markers = "PATH/TO/target-markers.fa"

# Global settings
#########################
# number of threads
threads = 16
# Amount of memory to allocate in GB
memory = 120
# TRUE to overwrite previous runs. FALSE the script will resume but will not delete anything.
overwrite = FALSE
# Hide verbose output for each function
quiet = FALSE

# Missing locus recovery settings
#########################
# TRUE = run expandMissingAssembly after the main assembly pipeline.
# Deprecated. Use the binned target assembly below instead.
expand.missing = FALSE
# Which read subdirectory within processed.reads to use for Phase 2 mapping.
# Should be paired (non-merged) reads. Options (use whichever is the last step run in workflow 1):
#   "decontaminated-reads" (default — recommended)
#   "cleaned-reads"
#   "error-corrected-reads"
# Do NOT use "pe-merged-reads" — HISAT2 expects paired input
mapping.reads = "decontaminated-reads"
# Reference to use for Phase 2 read mapping:
#   "contig"    = use the best assembled contig from other samples (default; closer match, better read recovery)
#   "reference" = use the original probe/bait sequences (useful when cross-sample contigs are absent or poor)
phase2.reference = "contig"
# TRUE = also attempt to recover loci absent from every sample's assembly,
#        mapping reads directly to the original reference sequences.
#        Can add substantial run time on large datasets.
recover.all.missing = FALSE
# Blast filters applied to the recovered contigs. These are usually less strict
# than the target contig filters above, because a recovered locus is expected to
# be shorter and more divergent.
expand.match.length = 100
expand.match.percent = 60
expand.match.coverage = 35

#Binned per-locus assembly settings
#########################
# TRUE = run assembleBinnedTargets after the main assembly pipeline.
# It assembles each target on its own, to recover targets the whole-library
# assembly lost and to extend the targets it found.
binned.assembly = FALSE
# Which read subdirectory within processed.reads to bin. Must be paired
# (non-merged) reads, because the mate of an anchored read supplies the flank.
binning.reads = "decontaminated-reads"
# The draft assembly directory that LAST searches for divergent contigs. A target
# assembles in the draft assembly whatever its divergence, but blastn cannot
# match the contig to a probe past about 25 percent divergence. LAST trains on
# the sample and finds those contigs, which then become the baits.
# Set to "" to skip the LAST search and use only the target contigs.
binned.draft.directory = "data-analysis/contigs/2_reduced-redundancy"
# Least part of the target length that a sample sequence must cover to be used as
# the bait, on a scale of 0 to 1. A shorter fragment recruits reads only across
# itself, so the reference is used instead.
binned.min.bait.coverage = 0.5
# TRUE = recruit reads with LAST for the targets that have no sequence in the
# sample, assemble them, and use the result as the bait. bwa cannot recruit a
# read that is 35 percent divergent from its bait, so without this step a
# divergent target that the draft assembly also lost cannot be recovered.
binned.rescue.missing = TRUE
# TRUE = after round 1, use LAST to recruit reads for the targets that no bin
# produced. bwa needs about 90 percent identity, so a divergent target recruits
# nothing. Costs one more pass over the reads, about 8 minutes per sample.
binned.rescue.failed.divergent = FALSE
# Which targets to bin:
#   "all"     = every target. Recovers missing loci and extends the ones present.
#   "missing" = only the targets absent from that sample. Much faster, and the
#               safe first run.
binned.locus.set = "all"
# What to map the reads against in the first round:
#   "hybrid"    = the contig of the sample where it has one, the reference for
#                 the rest. A contig of the sample recruits far more reads.
#   "reference" = always the target markers
binned.bait.source = "hybrid"
# Number of bait-and-assemble rounds. 1 adds about one insert length of flank to
# each side. Each later round adds about one more, at a growing risk of
# extension into a repeat.
binned.iterations = 1
# Minimum read pairs a bin needs before it is assembled
binned.min.pairs = 6
# Maximum read pairs kept per bin. 0 removes the limit.
binned.max.pairs = 3000
# Maximum base pairs a contig may add to each side of its bait in one round
binned.max.extension = 1000
# Maximum targets a new contig may match before it is treated as a repeat
binned.max.target.hits = 5
# What to do with a target that has more than one contig in the assembly:
#   "keep"    = leave those targets unchanged. One binned contig cannot stand
#               for two copies.
#   "longest" = treat them like any other target and keep the longest sequence.
binned.multi.copy = "keep"
# BLAST filters for the binned contigs. The coverage filter is lower than the
# target contig filter above, because a binned contig is already anchored to its
# own target.
binned.match.length = 50
binned.match.percent = 60
binned.match.coverage = 30
# The k-mer values for the per-bin SPAdes runs. Every value must be below the
# read length. Fewer values is faster.
binned.kmer.values = c(21, 33, 55, 77, 99)

#Contig curation settings
#########################
# TRUE = run curateTargetContigs after the steps above. It joins the fragments of
# one target that sit on separate contigs, and cuts apart a contig that spans
# more than one target. Run this before workflow 3, because the variant caller
# maps the reads back to these contigs.
# It runs twice: once after the target filter, so the steps below get joined
# single-locus baits, and once at the end, so workflow 3 maps reads to a clean set.
curate.contigs = TRUE
# Which program matches the target markers to the contigs. See target.search.method.
curate.search.method = "last"
# Match filters for the curation step. The coverage is summed over every hit of
# a target and tested after the fragments are joined, so a target split across
# two contigs is kept. The default is permissive on purpose: a later step can
# remove a short locus, but this step cannot recover one it dropped.
curate.match.length = 50
curate.match.percent = 60
curate.match.coverage = 30

#Assembly settings
#########################
#The selected k-mer values for spades
spades.kmer.values = c(33, 55, 77, 99, 127)
#Whether to use mismatch corrector (requires a lot of RAM and resources, recommended if possible)
spades.mismatch.corrector = TRUE
#TRUE runs spades in --isolate mode, for high coverage isolate data.
#Cannot be TRUE at the same time as spades.mismatch.corrector.
spades.isolate = FALSE
#Whether to save the error-corrected reads produced by SPAdes (ignored if clean.up.spades = TRUE)
save.corrected.reads = FALSE
#TRUE to delete the entire SPAdes working directory for each sample after assembly, keeping only the final .fa file
clean.up.spades = FALSE
# The similarity threshold for redundancy reduction. cd-hit-est needs a value of
# 0.8 or greater.
similarity = 0.95

#Target contig filtering settings
#########################
# A contig is kept when its blast hit to a target marker passes all three tests.
# Minimum blast alignment length in base pairs
target.match.length = 50
# Minimum blast percent identity, on a scale of 0 to 100
target.match.percent = 60
# Minimum percentage of the target marker length that the hit must cover
target.match.coverage = 30
# Which program matches the contigs to the target markers:
#   "last"  = LAST. It matches a contig up to about 35 percent divergent from
#             its target. Recommended.
#   "blast" = the previous blastn dc-megablast search. It loses a contig past
#             about 25 percent divergence.
target.search.method = "last"

#Program paths
#########################
### *** When installing the pipeline requirements via anaconda, only the path is needed to the conda bin directory
### Otherwise, if installed other ways, modify any of these to their path if R is not detecting system paths
conda.env = "/PATH/TO/miniconda3/envs/PhyloProcessR/bin"
cdhit.path    = conda.env
spades.path   = conda.env
megahit.path  = conda.env
cap3.path     = conda.env
blast.path    = conda.env
hisat2.path   = conda.env
bwa.path      = conda.env
last.path     = conda.env
samtools.path = conda.env
fastp.path    = conda.env
