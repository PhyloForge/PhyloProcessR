#' @title createContaminantDB
#'
#' @description Builds a local directory of contaminant reference genomes to be
#'   used by removeContamination(). Downloads genomes from GenBank via the
#'   biomartr package based on a user-supplied CSV, optionally adds the human
#'   genome, NCBI UniVec vector sequences, additional GenBank accessions, or a
#'   custom FASTA file.
#'
#' @param decontamination.list path to a CSV file with at least two columns:
#'   Genome (a short display name) and GenBank_Accession (the GenBank accession
#'   to download). Rows whose Genome matches "Human", "Homo", etc. can be
#'   excluded via include.human.
#'
#' @param output.directory path to the directory where contaminant genome files
#'   will be saved.
#'
#' @param include.human logical; if FALSE, rows associated with human genomes
#'   are removed from decontamination.list before downloading.
#'
#' @param include.univec logical; if TRUE, the NCBI UniVec vector/adaptor
#'   sequence database is downloaded and added to the contaminant directory.
#'
#' @param include.genbank character vector of additional GenBank organism names
#'   to download and include; NULL skips this step.
#'
#' @param include.fasta path to an existing FASTA file to copy directly into
#'   the contaminant directory; NULL skips this step.
#'
#' @param overwrite logical; if TRUE an existing output.directory is deleted
#'   and recreated. If FALSE and the directory already exists the function
#'   returns early with a message.
#'
#' @return invisibly; side effect is a populated output.directory containing
#'   compressed genome FASTA files ready for use as a decontamination reference.
#'
#' @export

createContaminantDB = function(decontamination.list = NULL,
                               output.directory = "contaminant-references",
                               include.human = TRUE,
                               include.univec = TRUE,
                               include.genbank = NULL,
                               include.fasta = NULL,
                               overwrite = FALSE) {

  # #Debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # output.directory = "contaminant-references"
  # decontamination.list = "decontamination_database.csv"
  # overwrite = FALSE
  # quiet = TRUE
  # include.human = TRUE
  # include.univec = TRUE
  # include.genbank = NULL
  # include.fasta = NULL
  # # map.match = 0.99

  #Quick checks
  if (is.null(decontamination.list) == TRUE){ stop("Please provide list of sequences and genbank numbers.") }
  if (file.exists(decontamination.list) == F){ stop("Input list not found.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    } else { return("Contaminant database already exists, skipping.") }
  }#end else

  sample.data = read.csv(file = decontamination.list)

  if (include.human == FALSE){ sample.data = sample.data[!sample.data$Genome %in% c("Human", "Homo", "Sapien", "Homo_Sapien", "Homo Sapien"),] }

  if (nrow(sample.data) == 0){ stop("no samples remain to analyze.") }

  if (is.null(include.fasta) != TRUE){
    if (file.exists(include.fasta) == F){ stop("include.fasta file not found.") }
    system(paste0("cp ", include.fasta, " ", output.directory, "/manually-included-data.fa"))
  }#end if

  # Downloads a genome by accession using the NCBI Datasets API (GCA/GCF) or
  # Entrez efetch (nucleotide accessions like NC_*), saves as a gzipped FASTA.
  .download_accession = function(accession, out.path) {
    is.assembly = grepl("^GC[AF]_", accession)
    if (is.assembly) {
      zip.file = tempfile(fileext = ".zip")
      url = paste0(
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/",
        accession,
        "/download?include_annotation_type=GENOME_FASTA"
      )
      download.file(url, destfile = zip.file, quiet = TRUE, mode = "wb")
      tmp.dir = tempfile()
      dir.create(tmp.dir)
      utils::unzip(zip.file, exdir = tmp.dir)
      fna.file = list.files(tmp.dir, pattern = "\\.fna$", recursive = TRUE, full.names = TRUE)
      if (length(fna.file) == 0) stop(paste("No .fna found in NCBI zip for", accession))
      system(paste0("gzip -c ", fna.file[1], " > ", shQuote(out.path)))
      unlink(zip.file); unlink(tmp.dir, recursive = TRUE)
    } else {
      # Nucleotide accession (NC_*, AY_*, etc.) via Entrez efetch
      url = paste0(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        "?db=nucleotide&id=", accession, "&rettype=fasta&retmode=text"
      )
      fa.file = tempfile(fileext = ".fa")
      download.file(url, destfile = fa.file, quiet = TRUE)
      system(paste0("gzip -c ", shQuote(fa.file), " > ", shQuote(out.path)))
      unlink(fa.file)
    }
  }

  if (is.null(include.genbank) != TRUE) {
    for (i in 1:length(include.genbank)) {
      out.path = file.path(output.directory, paste0("include-genbank_", include.genbank[i], ".fna.gz"))
      message("Downloading include.genbank: ", include.genbank[i])
      tryCatch(
        .download_accession(include.genbank[i], out.path),
        error = function(e) warning("Could not download ", include.genbank[i], ": ", conditionMessage(e))
      )
    }
  }

  if (include.univec == TRUE) {
    download.file(
      "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec",
      destfile = file.path(output.directory, "UniVec.fa"),
      quiet = TRUE
    )
  }
  sample.data = sample.data[!sample.data$Genome %in% c("univec", "UniVec", "UNIVEC", "Univec"), ]

  for (i in 1:nrow(sample.data)) {
    accession = trimws(sample.data$GenBank_Accession[i])
    # Skip rows whose accession column contains a URL (handled elsewhere)
    if (grepl("^https?://", accession)) next
    out.path = file.path(output.directory,
                         paste0(sample.data$Genome[i], "-", accession, ".fna.gz"))
    message("Downloading ", sample.data$Genome[i], " (", accession, ")")
    tryCatch(
      .download_accession(accession, out.path),
      error = function(e) warning("Could not download ", accession, ": ", conditionMessage(e))
    )
  }

  file.list = list.files(output.directory, full.names = TRUE)
  del.files = file.list[grep("\\.fna\\.gz$|\\.fa$", file.list, invert = TRUE)]
  if (length(del.files) > 0) {
    unlink(del.files, recursive = TRUE)
  }

  final.files = list.files(output.directory)

  if (length(final.files) != 0){ print("Creation of contamination database was successful.") }

} #end function

