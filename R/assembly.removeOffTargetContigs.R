#' @title removeOffTargetContigs
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param input.reads path to a folder of raw reads in fastq format.
#'
#' @param target.markers a csv file with a "File" and "Sample" columns, where "File" is the file name and "Sample" is the desired renamed file
#'
#' @param output.name the new directory to save the adaptor trimmed sequences
#'
#' @param mapper "Sample" to run on a single sample or "Directory" to run on a directory of samples
#'
#' @param min.iterations system path to fastp in case it can't be found
#'
#' @param max.iterations system path to fastp in case it can't be found
#'
#' @param min.length system path to fastp in case it can't be found
#'
#' @param max.length system path to fastp in case it can't be found
#'
#' @param min.ref.id system path to fastp in case it can't be found
#'
#' @param spades.path system path to fastp in case it can't be found
#'
#' @param bbmap.path system path to fastp in case it can't be found
#'
#' @param cap3.path system path to fastp in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
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

#Iteratively assembles to reference
removeOffTargetContigs = function(assembly.directory = NULL,
                                  target.markers = NULL,
                                  output.directory = "target-contigs",
                                  blast.path = NULL,
                                  memory = 1,
                                  threads = 1,
                                  overwrite = FALSE,
                                  quiet = TRUE) {

  #Debug
  # library(PhyloProcessR)
  # setwd("/Volumes/LaCie/Mantellidae")
  # assembly.directory = "data-analysis/contigs/reduced-redundancy"
  # output.directory = "data-analysis/contigs/target-contigs"
  # target.markers = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
  # blast.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  # quiet = TRUE
  # overwrite = FALSE
  # threads = 8
  # memory = 20

  if (is.null(blast.path) == FALSE) {
    b.string <- unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    blast.path <- ""
  }

  # Quick checks
  if (is.null(assembly.directory) == TRUE) {
    stop("Please provide input reads.")
  }
  if (is.null(target.markers) == TRUE) {
    stop("Please provide a reference.")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else {
    dir.create(output.directory)
  }

  file.names = list.files(assembly.directory)
 
  #############################
  ## Target matching loop start
  #############################

  #headers
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  for (i in seq_along(file.names)) {
    # Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])
    species.dir = paste0(output.directory, "/", sample)

    # Creates species directory if none exists
    if (file.exists(species.dir) == FALSE) {
      dir.create(species.dir)
    }

    # Make blast database for the probe loci
    system(paste0(
      blast.path, "makeblastdb -in ", target.markers,
      " -parse_seqids -dbtype nucl -out ", species.dir, "/nucl-blast_db"
    ), ignore.stdout = quiet)

    # Matches samples to loci
    system(paste0(
      blast.path, "blastn -task dc-megablast -db ", species.dir, "/nucl-blast_db -evalue 0.001",
      " -query ", assembly.directory, "/", file.names[i], " -out ", species.dir, "/target-blast-match.txt",
      " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
      " -num_threads ", threads
    ))

    # Loads in match data
    match.data = data.table::fread(paste0(species.dir, "/target-blast-match.txt"),
      sep = "\t", header = FALSE, stringsAsFactors = FALSE
    )
    data.table::setnames(match.data, headers)

    # Matches need to be greater than 12
    filt.data = match.data[match.data$matches > 60, ]
    # Percent identitiy must match 50% or greater
    filt.data = filt.data[filt.data$pident >= 0.6, ]

    # Make sure the hit is greater than 50% of the target.markers length
    filt.data = filt.data[filt.data$matches >= ((30 / 100) * filt.data$tLen), ]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      next
    }

    # Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    dup.contigs = unique(filt.data[duplicated(filt.data$tName) == TRUE, ]$tName)
    good.data = filt.data[!filt.data$tName %in% dup.contigs, ]

    dedup.data = c()
    for (j in seq_along(dup.contigs)){
        #remove duplicates
        temp.data = filt.data[filt.data$tName == dup.contigs[j], ]
        temp.data = temp.data[duplicated(temp.data$qName) == FALSE, ]
        dedup.data = rbind(dedup.data, temp.data)
    }#end j loop

    new.data = rbind(good.data, dedup.data)

    #########################################################################
    # Part B: Multiple sample contigs (tName) matching to one target (qName)
    #########################################################################
    # Pulls out
    og.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", file.names[i]),
      format = "fasta"
    )

    save.contigs = c()
    for (j in seq_len(nrow(new.data))) {
      temp.contigs = og.contigs[names(og.contigs) %in% new.data$qName[j]]
      names(temp.contigs) = new.data$tName[j]
      save.contigs = append(save.contigs, temp.contigs)
    } # end j loop

    # Renames duplicates
    names(save.contigs) = make.unique(names(save.contigs), sep = "_")

    # Finds probes that match to two or more contigs
    final.loci = as.list(as.character(save.contigs))
    writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", sample, ".fa"), nbchar = 1000000, as.string = TRUE
    )

    system(paste0("rm -r ", species.dir))
  } # end iterations if

  ##########################
}#end function