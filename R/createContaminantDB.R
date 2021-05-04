#' @title createContaminantDB
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param decontamination.list path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param include.human the new directory to save the adaptor trimmed sequences
#'
#' @param include.mouse directory of genomes contaminants to scan samples
#'
#' @param overwrite system path to samtools in case it can't be found
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export

createContaminantDB = function(decontamination.list = NULL,
                               output.directory = "contaminant-references",
                               include.human = TRUE,
                               include.univec = TRUE,
                               overwrite = FALSE) {

  #Debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # read.directory = "read-processing/adaptor-removed-reads"
  # output.directory = "read-processing/decontaminated-reads"
  # decontamination.path = "/Users/chutter/Dropbox/Research/0_Github/Contamination_Genomes"
  # samtools.path = "/Users/chutter/miniconda3/bin/samtools"
  # bwa.path = "/usr/local/bin/bwa"
  # mode = "directory"
  # threads = 4
  # mem = 8
  # resume = TRUE
  # overwrite = FALSE
  # quiet = TRUE
  # map.match = 0.99

  #Quick checks
  if (is.null(decontamination.list) == TRUE){ stop("Please provide list of sequences and genbank numbers.") }
  if (file.exists(decontamination.list) == F){ stop("Input list not found.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  sample.data = read.csv(file = decontamination.list)

  if (include.human == FALSE){ sample.data = sample.data[!sample.data$Genome %in% c("Human", "Homo", "Sapien", "Homo_Sapien", "Homo Sapien"),] }

  if (nrow(sample.data) == 0){ stop("no samples remain to analyze.") }

  if (include.univec == TRUE){
   system(paste0("wget -c https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec -O ", output.directory, "/UniVec.fa"))
  }
  sample.data = sample.data[!sample.data$Genome %in% c("univec", "UniVec", "UNIVEC", "Univec"),]


  for (i in 1:nrow(sample.data)){
    biomartr::getGenome(db = "genbank",
                        organism = sample.data$GenBank_Accession[i],
                        path = output.directory,
                        reference = FALSE)

    new.name = paste0(sample.data$Genome[i], "-", sample.data$GenBank_Accession[i], ".fna.gz")

    file.list = list.files(output.directory)
    sample.files = file.list[grep(sample.data$GenBank_Accession[i], file.list)]
    old.name = sample.files[grep(".fna.gz|.fa$", sample.files)]

    system(paste0("mv ", output.directory, "/", old.name, " ", output.directory, "/", new.name))

  }#end i loop

  file.list = list.files(output.directory)
  del.files = file.list[grep(".fna.gz|.fa$", file.list, invert = T)]
  system(paste0("rm ", paste0(output.directory, "/", del.files, collapse = " ")))

  final.files = list.files(output.directory)

  if (length(final.files) != 0){ print("Creation of contamination database was successful.") }

} #end function

