#' @title removeOffTargetContigs
#'
#' @description Filters each sample's assembly to retain only contigs that have
#'   a BLAST match to the provided target markers. One BLAST database is built
#'   from the target markers and each assembly is queried against it
#'   (dc-megablast). Hits are filtered by alignment length, percent identity, and
#'   coverage of the target. Surviving contigs are renamed to match their target
#'   and saved as a single FASTA per sample in the output directory. When several
#'   contigs match the same target, each contig is kept once, with its best hit,
#'   and the duplicate names are made unique with a numeric suffix. Samples that
#'   already have a FASTA file in the output directory are skipped unless
#'   \code{overwrite = TRUE}.
#'
#' @param assembly.directory path to a directory of assembly FASTA files (.fa),
#'   one per sample.
#'
#' @param target.markers path to a FASTA file of target marker sequences used
#'   as the BLAST reference database.
#'
#' @param output.directory path to the directory where filtered per-sample FASTA
#'   files will be written. Default: \code{"target-contigs"}.
#'
#' @param search.method which program matches the contigs to the target markers.
#'   \code{"last"} (default) uses LAST, which finds a contig that is up to about
#'   35 percent divergent from its target. \code{"blast"} uses the previous
#'   \code{blastn dc-megablast} search, which loses a contig past about 25
#'   percent divergence. Numbers in HANDOFF.md.
#'
#' @param last.path path to the directory containing \code{lastdb} and
#'   \code{lastal}. Only needed when \code{search.method = "last"}. If
#'   \code{NULL} the programs must be on the system PATH. Default: \code{NULL}.
#'
#' @param blast.path path to the directory containing \code{makeblastdb} and
#'   \code{blastn}. If \code{NULL} expected on the system PATH. Default:
#'   \code{NULL}.
#'
#' @param min.match.length minimum BLAST alignment length (bp) a hit must exceed
#'   to be accepted. Default: \code{60}.
#'
#' @param min.match.percent minimum BLAST percent identity (0-100) required to
#'   accept a hit. Default: \code{60}.
#'
#' @param min.match.coverage minimum percentage of the target length that a
#'   BLAST hit must cover to be accepted. Default: \code{30}.
#'
#' @param memory not currently used; reserved for future use. Default: \code{1}.
#'
#' @param threads number of CPU threads passed to \code{blastn}. Default:
#'   \code{1}.
#'
#' @param overwrite logical; if \code{TRUE} the output directory is deleted and
#'   recreated and all samples are rerun. Default: \code{FALSE}.
#'
#' @param quiet logical; if \code{TRUE} \code{makeblastdb} screen output is
#'   suppressed. Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Writes one filtered FASTA file per sample
#'   to \code{output.directory}, with contigs renamed to match their target
#'   marker names.
#'
#' @export

#Iteratively assembles to reference
removeOffTargetContigs = function(assembly.directory = NULL,
                                  target.markers = NULL,
                                  output.directory = "target-contigs",
                                  blast.path = NULL,
                                  last.path = NULL,
                                  search.method = c("last", "blast"),
                                  min.match.length = 60,
                                  min.match.percent = 60,
                                  min.match.coverage = 30,
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

  search.method = match.arg(search.method)
  last.path = .programPrefix(last.path)

  # Quick checks
  if (is.null(assembly.directory) == TRUE) {
    stop("Please provide input reads.")
  }
  if (is.null(target.markers) == TRUE) {
    stop("Please provide a reference.")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) {
      unlink(output.directory, recursive = TRUE)
      dir.create(output.directory, recursive = TRUE)
    }
  } else {
    dir.create(output.directory, recursive = TRUE)
  }

  # Only FASTA files are samples. Hidden files and other output are ignored.
  fasta.pattern = "\\.fa$|\\.fas$|\\.fasta$|\\.fna$"
  file.names = list.files(assembly.directory, pattern = fasta.pattern)

  # Resume: skip samples already written to the output directory
  if (overwrite == FALSE) {
    done.names = gsub(fasta.pattern, "", list.files(output.directory, pattern = fasta.pattern))
    file.names = file.names[!gsub(fasta.pattern, "", file.names) %in% done.names]
  }

  if (length(file.names) == 0) { return(invisible(NULL)) }

  # Build the search database once from the target markers. Every sample uses it.
  blast.db.dir = file.path(output.directory, "blast_db")
  dir.create(blast.db.dir, showWarnings = FALSE, recursive = TRUE)

  if (search.method == "last") {
    .lastBuildDB(reference.file = target.markers,
                 db.prefix = paste0(blast.db.dir, "/last_db"),
                 lastdb.command = paste0(last.path, "lastdb"),
                 threads = threads,
                 quiet = quiet)
  } else {
    system(paste0(
      blast.path, "makeblastdb -in ", shQuote(target.markers),
      " -parse_seqids -dbtype nucl -out ", shQuote(paste0(blast.db.dir, "/nucl-blast_db"))
    ), ignore.stdout = quiet)
  }

  #############################
  ## Target matching loop start
  #############################

  #headers
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  for (i in seq_along(file.names)) {
    # Sets up working directories for each species
    sample = gsub(pattern = fasta.pattern, replacement = "", x = file.names[i])
    species.dir = paste0(output.directory, "/", sample)
    blast.file = paste0(species.dir, "/target-blast-match.txt")

    # Creates species directory if none exists
    if (dir.exists(species.dir) == FALSE) {
      dir.create(species.dir, recursive = TRUE)
    }

    # Matches samples to loci
    if (search.method == "last") {
      .lastSearch(query.file = paste0(assembly.directory, "/", file.names[i]),
                  db.prefix = paste0(blast.db.dir, "/last_db"),
                  out.file = blast.file,
                  lastal.command = paste0(last.path, "lastal"),
                  threads = threads,
                  quiet = quiet)
    } else {
      system(paste0(
        blast.path, "blastn -task dc-megablast -db ", shQuote(paste0(blast.db.dir, "/nucl-blast_db")),
        " -evalue 0.001",
        " -query ", shQuote(paste0(assembly.directory, "/", file.names[i])),
        " -out ", shQuote(blast.file),
        " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
        " -num_threads ", threads
      ))
    }

    # An empty file means no hits. fread returns no columns for an empty file, so
    # the file is checked before the names are set.
    if (file.exists(blast.file) == FALSE || file.size(blast.file) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      unlink(species.dir, recursive = TRUE)
      next
    }

    # Loads in match data
    match.data = data.table::fread(blast.file,
      sep = "\t", header = FALSE, stringsAsFactors = FALSE
    )
    data.table::setnames(match.data, headers)

    # Alignment length must be greater than min.match.length
    filt.data = match.data[match.data$matches > min.match.length, ]
    # Percent identity must be at least min.match.percent
    filt.data = filt.data[filt.data$pident >= min.match.percent, ]

    # Make sure the hit covers at least min.match.coverage of the target length
    filt.data = filt.data[filt.data$matches >= ((min.match.coverage / 100) * filt.data$tLen), ]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      unlink(species.dir, recursive = TRUE)
      next
    }

    # Sorting: contig name, target name, bitscore higher first, identity, evalue.
    # Bitscore comes before identity so a long strong hit beats a short exact one.
    data.table::setorder(filt.data, qName, tName, -bitscore, -pident, evalue)

    dup.targets = unique(filt.data[duplicated(filt.data$tName) == TRUE, ]$tName)
    good.data = filt.data[!filt.data$tName %in% dup.targets, ]

    dedup.data = c()
    for (j in seq_along(dup.targets)){
        #keeps the best hit for each contig that matches this target
        temp.data = filt.data[filt.data$tName == dup.targets[j], ]
        temp.data = temp.data[duplicated(temp.data$qName) == FALSE, ]
        dedup.data = rbind(dedup.data, temp.data)
    }#end j loop

    new.data = rbind(good.data, dedup.data)

    #########################################################################
    # Part B: Multiple sample contigs (qName) matching to one target (tName)
    #########################################################################
    og.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", file.names[i]),
      format = "fasta"
    )

    # BLAST reports the first word of the header, so the contig names are cut at
    # the first space before they are matched.
    names(og.contigs) = gsub(" .*", "", names(og.contigs))

    # One vectorised lookup instead of growing the set one contig at a time
    keep.index = match(new.data$qName, names(og.contigs))
    new.data = new.data[!is.na(keep.index), ]
    keep.index = keep.index[!is.na(keep.index)]

    if (length(keep.index) == 0) {
      print(paste0(sample, " had no matching contigs. Skipping"))
      unlink(species.dir, recursive = TRUE)
      next
    }

    save.contigs = og.contigs[keep.index]

    # Renames each contig to its target and makes duplicate names unique
    names(save.contigs) = make.unique(new.data$tName, sep = "_")

    final.loci = as.list(as.character(save.contigs))
    writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", sample, ".fa"), nbchar = 1000000, as.string = TRUE
    )

    unlink(species.dir, recursive = TRUE)
  } # end iterations if

  unlink(blast.db.dir, recursive = TRUE)

}#end function