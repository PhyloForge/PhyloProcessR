

#Manual installation
install.packages(c("devtools", "ape", "data.table", "doparallel", "foreach", "iterators", "rdrop2", "seqinr", "stringr"),
                 dependencies = T, repos = "http://cran.us.r-project.org", upgrade = "never", force = TRUE)


#devtools::install_version("ape", version = "5.4.1", repos = "http://cran.us.r-project.org",
#                          upgrade = "never", force = TRUE)
#devtools::install_version("data.table", version = "1.13", repos = "http://cran.us.r-project.org")
#devtools::install_version("doparallel", version = "1.0.16", repos = "http://cran.us.r-project.org")
#devtools::install_version("foreach", version = "1.5.1", repos = "http://cran.us.r-project.org")
#devtools::install_version("iterators", version = "1.0.13", repos = "http://cran.us.r-project.org")
#devtools::install_version("rdrop2", version = "0.8.2.1", repos = "http://cran.us.r-project.org")
#devtools::install_version("seqinr", version = "4.2", repos = "http://cran.us.r-project.org")
#devtools::install_version("stringr", version = "1.4.0", repos = "http://cran.us.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install(c("biostrings", "genomicranges", "rsamtools"))

#Installs updated package version
devtools::install_github("chutter/PhyloCap", upgrade = "never")
library(PhyloCap)
library(foreach)

threads = 4
memory = 8

#################################################
## Configuration file
##################

#Directories
setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
read.dir = "raw-reads"
decontamination.path = "/Users/chutter/Dropbox/Research/0_Github/Contamination_Genomes"
file.rename = "/Users/chutter/Dropbox/Research/0_Github/Test-dataset/file_rename.csv"
target.file = "/Users/chutter/Dropbox/Research/0_Github/PhyloCap/setup-configuration_files/Ranoidea_All-Markers_Apr21-2019.fa"
dataset.name = "Test"

#Pre-processing settings
decontamination = TRUE
decontamination.match = 0.99 #
spades.kmer.values = c(21,33,55,77,99,127)
spades.mismatch.corrector = TRUE


#Target matching and alignment settings
min.percent.id = 0.5
min.match.length = 40
min.match.coverage = 0.5
min.taxa.alignment = 4 #minimum taxa for an alignment


#Program paths
#conda.bin.path = "/Users/chutter/conda/PhyloCap/bin"
fastp.path = "/Users/chutter/conda/PhyloCap/bin"
samtools.path = "/Users/chutter/conda/PhyloCap/bin"
bwa.path = "/Users/chutter/conda/PhyloCap/bin"
spades.path = "/Users/chutter/conda/PhyloCap/bin"
bbmap.path = "/Users/chutter/conda/PhyloCap/bin"
blast.path = "/Users/chutter/conda/PhyloCap/bin"
mafft.path = "/Users/chutter/conda/PhyloCap/bin"
iqtree.path = "/Users/chutter/conda/PhyloCap/bin"

### Example usage
#setwd("/home/c111h652/scratch/Shrew_UCE")
#read.dir = "/home/c111h652/scratch/Shrew_UCE"
#decontamination.path = "/home/c111h652/scratch/Contamination_Genomes"

#fastp.path = "/panfs/pfs.local/work/bi/c111h652/conda/phylocap/bin/fastp"
#samtools.path = "/panfs/pfs.local/work/bi/c111h652/conda/phylocap/bin/samtools"
#bwa.path = "/panfs/pfs.local/work/bi/c111h652/conda/phylocap/bin/bwa"
#spades.path = "/panfs/pfs.local/work/bi/c111h652/conda/phylocap/bin/spades.py"
#bbmap.path = "/panfs/pfs.local/work/bi/c111h652/conda/phylocap/bin"
#iqtree.path = "/home/c111h652/programs"

#################################################
## Step 2: Clean out contamination
##################

dir.create("processed-reads")

organizeReads(read.directory = read.dir,
              output.dir = "processed-reads/organized-reads",
              rename.file = file.rename,
              overwrite = FALSE)

removeAdaptors(input.reads = "processed-reads/organized-reads",
               output.directory = "processed-reads/adaptor-removed-reads",
               fastp.path = fastp.path,
               threads = threads,
               mem = memory,
               resume = FALSE,
               overwrite = TRUE,
               quiet = TRUE)

## remove external contamination
removeContamination(input.reads = "processed-reads/adaptor-removed-reads",
                    output.directory = "processed-reads/decontaminated-reads",
                    decontamination.path = decontamination.path,
                    map.match = decontamination.match,
                    samtools.path = samtools.path,
                    bwa.path = bwa.path,
                    threads = threads,
                    mem = memory,
                    resume = FALSE,
                    overwrite = TRUE,
                    overwrite.reference = FALSE,
                    quiet = FALSE)

#merge paired end reads
mergePairedEndReads(input.reads = "processed-reads/decontaminated-reads",
                    output.directory = "processed-reads/pe-merged-reads",
                    fastp.path = fastp.path,
                    threads = threads,
                    mem = memory,
                    resume = FALSE,
                    overwrite = TRUE,
                    quiet = TRUE)

#Assembles merged paired end reads with spades
assembleSpades(input.reads = "processed-reads/pe-merged-reads",
               output.directory = "processed-reads/spades-assembly",
               assembly.directory = "draft-assemblies",
               mismatch.corrector = spades.mismatch.corrector,
               kmer.values = spades.kmer.values,
               threads = threads,
               memory = memory,
               overwrite = FALSE,
               resume = TRUE,
               save.corrected.reads = FALSE,
               quiet = TRUE,
               spades.path = spades.path)

#match targets
matchTargets(assembly.directory = "draft-assemblies",
             target.file = target.file,
             alignment.contig.name = dataset.name,
             output.directory = "match-targets",
             min.percent.id = 0.5,
             min.match.length = 40,
             min.match.coverage = 0.5,
             threads = threads,
             memory = memory,
             trim.target = FALSE,
             overwrite = TRUE,
             resume = FALSE,
             quiet = TRUE,
             blast.path = blast.path,
             bbmap.path = bbmap.path)

#align targets
alignTargets(targets.to.align = paste0(dataset.name, "_match-targets_to-align.fa"),
             output.directory = "alignments",
             min.taxa = min.taxa.alignment,
             subset.start = 0,
             subset.end = 1,
             threads = threads,
             memory = memory,
             overwrite = FALSE,
             resume = TRUE,
             quiet = TRUE,
             mafft.path = mafft.path)



#Fix the installs for this
batchTrimAlignments(alignment.dir = paste0("alignments"),
                    alignment.format = "phylip",
                    output.dir = paste0("alignments-trimmed"),
                    output.format = "phylip",
                    overwrite = FALSE,
                    resume = TRUE,
                    TAPER = TRUE,
                    TAPER.path = "/home/c111h652/programs/correction.jl",
                    julia.path = "julia",
                    TrimAl = TRUE,
                    HmmCleaner = FALSE,
                    trim.column = TRUE,
                    convert.ambiguous.sites = TRUE,
                    alignment.assess = TRUE,
                    trim.external = TRUE,
                    trim.coverage = TRUE,
                    min.coverage.percent = 35,
                    min.external.percent = 50,
                    min.column.gap.percent = 50,
                    min.align.length = 100,
                    min.taxa.count = 5,
                    min.gap.percent = 50,
                    min.sample.bp = 60,
                    threads = threads,
                    memory = memory)


estimateGeneTrees(alignment.directory = "alignments-trimmed",
                  output.directory = "gene-trees",
                  min.taxa = 4,
                  subset.start = 0,
                  subset.end = 1,
                  threads = threads,
                  memory = memory,
                  overwrite = FALSE,
                  resume = TRUE,
                  quiet = TRUE,
                  cleanup.files = TRUE,
                  iqtree.path = iqtree.path)


