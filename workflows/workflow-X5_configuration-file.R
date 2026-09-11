#################################################
## Configuration file for workflow X5: paralog analysis
#################################################

# Package version
#########################
# TRUE installs the latest development version from GitHub before the run.
install.latest.github = FALSE

# Directories and input files
#########################
working.directory = "/PATH/TO/PROJECT/DIRECTORY"
alignment.directory = "data-analysis/alignments/untrimmed_all-markers"
alignment.format = "phylip"
paralog.directory = "data-analysis/contigs/9_paralog-contigs"
primary.directory = "data-analysis/contigs/8_annotated-contigs"
candidate.directory = "logs/sample_logs"
target.file = "/PATH/TO/CAPTURE-TARGETS.fa"
output.directory = "data-analysis/paralog-analysis"

# NULL processes the full union of base-alignment and saved-paralog targets.
# Supply exact target IDs only for a labeled pilot or development subset.
target.names = NULL

# Optional curated references
#########################
# Supply both files or keep both values NULL.
reference.file = NULL
reference.table = NULL

# Stage controls
#########################
# A disabled stage requires compatible completed outputs from an earlier run.
collect.copies = TRUE
align.copies = TRUE
trim.copies = TRUE
infer.trees = TRUE
separate.copies = TRUE

# Global settings
#########################
threads = 8
memory = 40
overwrite = FALSE
quiet = TRUE

# Alignment settings
#########################
alignment.algorithm = "localpair"

# Copy-aware trimming settings
#########################
run.TrimAl = TRUE
min.external.percent = 50
min.column.gap.percent = 50
min.coverage.percent = 35
min.coverage.bp = 60
min.alignment.length = 100
min.taxa.alignment = 4
max.alignment.gap.percent = 50

# Tree and split settings
#########################
tree.model = "MFP"
bootstrap.replicates = 1000
tree.seed = 12345
min.branch.support = 95
min.shared.samples = 3
min.shared.sample.fraction = 0.5
min.split.branch.length = 0.05
min.split.branch.ratio = 10
# Use "retain" only to retain assessed single-copy targets with review flags.
review.action = "exclude"

# Optional downstream datasets
#########################
build.downstream.datasets = FALSE
feature.gene.names = NULL
minimum.exons = 2
# Optional TSV columns: Output_marker, Gene_copy_id, Unlinked_eligible.
verified.copy.mapping = NULL

# Program paths
#########################
conda.env = "PATH/TO/miniconda3/envs/PhyloProcessR/bin"
mafft.path = conda.env
trimAl.path = conda.env
iqtree.path = conda.env
iqtree.executable = "iqtree2"

#### End configuration
