#################################################
## Configuration file for workflow 6: legacy integration
#################################################

# Directories and input files
#########################
# *** Full paths should be used whenever possible
# The main working directory
working.directory = "/PATH/TO/PROJECT/DIRECTORY"
# The sequence capture target marker file used to BLAST-match each legacy locus
target.file = "/PATH/TO/marker-seqs.fa"
# Feature/gene name metadata file: column 1 "Marker", column 2 "Gene"
feature.gene.names = "data-analysis/gene_metadata.txt"

# Input alignment directory (default: untrimmed all-markers from workflow 4/5)
alignment.directory = "data-analysis/alignments/untrimmed_all-markers"
alignment.format = "phylip"

# Nexus conversion (optional — run before legacy integration)
#########################
# TRUE = split a concatenated NEXUS file into per-locus phylip files before integration.
#        The converted files become the legacy alignments used in integrateLegacy.
convert.nexus = FALSE
# Full path to the concatenated NEXUS file containing a BEGIN SETS / charset block
nexus.file = NULL
# Directory to write the per-locus phylip files produced by the conversion.
# When convert.nexus = TRUE this is used automatically as legacy.directory.
nexus.output.directory = "data-analysis/legacy-alignments"
# Maximum percent of missing/gap characters allowed per sample per locus.
# Samples exceeding this are dropped from that locus. Default 100 = keep all.
max.missing.percent = 100

# Legacy (Sanger/GenBank) alignment directory
# If convert.nexus = TRUE this is set automatically to nexus.output.directory
legacy.directory = "/PATH/TO/legacy-alignments"
legacy.format = "phylip"

# Global settings
#########################
# Number of threads
threads = 8
# Amount of memory to allocate in GB
memory = 40
# Whether to overwrite previous runs
overwrite = FALSE
# Hide verbose output for each function
quiet = TRUE

# Legacy integration settings
#########################
# TRUE = merge sequences from the same sample found in both datasets into one sequence
combine.same.sample = TRUE
# How to match sample names between legacy and capture alignments:
#   "exact"   = names must be identical
#   "species" = strip trailing specimen/voucher ID before matching;
#               merged sequence is named to species only (e.g. Genus_species)
name.match = "exact"
# TRUE = also include legacy loci absent from the capture dataset as stand-alone alignments
include.uncaptured.legacy = FALSE
# TRUE = write a second output directory containing all capture alignments updated
#        with legacy-integrated versions where available (recommended)
include.all.together = TRUE

# Mitochondrial loci settings (optional)
#########################
# TRUE = also integrate mitochondrial legacy loci that are absent from the nuclear
#        target marker file. A second BLAST database is built from consensus sequences
#        of the mitochondrial capture alignments in mito.alignment.directory.
#        Recommended: generate mito capture alignments with MitoTrawlR
#        (https://github.com/chutter/MitoTrawlR).
include.mitochondrial = FALSE
# Directory containing the mitochondrial capture alignment files.
# Only used when include.mitochondrial = TRUE.
mito.alignment.directory = "/PATH/TO/mito-alignments"
# Format of the mitochondrial alignment files ("phylip" or "fasta")
mito.alignment.format = "phylip"

# Gene concatenation settings
#########################
# TRUE = concatenate exons from the same gene using feature.gene.names metadata
concatenate.genes = TRUE
# Minimum number of exons needed to produce a concatenated gene alignment
minimum.exons = 2
# TRUE  = include legacy-only loci in the gene concatenation step
# FALSE = concatenate only capture-derived exons; legacy loci remain as separate alignments
concatenate.legacy.genes = TRUE
# TRUE = gather all unlinked alignments (concatenated genes + single-exon loci) for analysis
gather.unlinked = TRUE

# Trimming alignment settings
#########################
# TRUE = run alignment trimming (superTrimmer)
trim.alignments = TRUE
# Minimum number of taxa to keep an alignment
min.taxa.alignment = 4
# Minimum alignment length in basepairs to keep an alignment
min.alignment.length = 100
# Run TrimAl to remove highly variable or misaligned columns
run.TrimAl = TRUE
# TRUE = trim columns below a minimum gap threshold
trim.column = TRUE
# Minimum percent of bases that must be present to keep a column
min.column.gap.percent = 50
# Resolve ambiguous IUPAC sites to a random constituent base
convert.ambiguous.sites = FALSE
# TRUE = trim poorly covered edges of each alignment
trim.external = TRUE
# Minimum percent of bases required to keep an edge column
min.external.percent = 50
# TRUE = remove samples below a minimum coverage threshold
trim.coverage = TRUE
# Minimum percent of bases a sample must have across the alignment
min.coverage.percent = 35
# Minimum number of bases a sample must have across the alignment
min.coverage.bp = 60

# Program paths
#########################
### *** When installing via conda, only the path to the conda bin directory is needed.
### Otherwise modify any of these to the path where the program is found.
conda.env = "PATH/TO/miniconda3/envs/PhyloProcessR/bin"
blast.path = conda.env
mafft.path = conda.env
trimAl.path = conda.env

#### End configuration
