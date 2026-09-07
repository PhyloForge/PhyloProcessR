#' @title removeContamination
#'
#' @description Removes reads that map to a set of contaminant genomes (e.g.
#'   human, mouse, or vector sequences) using BWA-MEM and samtools. A combined
#'   BWA index is built from all FASTA files in decontamination.path and the
#'   reads are mapped against it. A read pair is a contaminant when either mate
#'   aligns at or above the map.match identity threshold. Contaminant pairs are
#'   removed from the read set and counted per contaminant genome. Every other
#'   pair is kept, including a pair that aligns below the threshold.
#'
#' @param input.reads path to a directory of adapter-trimmed paired-end reads
#'   in fastq.gz format, organised in per-sample sub-directories.
#'
#' @param output.directory path to the directory where decontaminated reads
#'   will be saved (one sub-directory per sample).
#'
#' @param decontamination.path path to a directory of contaminant reference
#'   genome FASTA files (e.g. produced by createContaminantDB()).
#'
#' @param map.match numeric between 0 and 1; the minimum alignment identity that
#'   makes a read a contaminant. 0.90 means 90 percent identity. Identity is
#'   calculated across the aligned CIGAR operations M, I, D, =, and X, excluding
#'   soft-clipped bases. The NM tag counts mismatches and inserted/deleted bases.
#'   The value sets both which pairs are removed and which pairs are counted, so
#'   the read set and the contamination report always agree.
#'
#' @param samtools.path system path to the directory that contains the samtools
#'   executable, or the full path to the executable; NULL searches the system
#'   PATH.
#'
#' @param bwa.path system path to the directory that contains the bwa
#'   executable, or the full path to the executable; NULL searches the system
#'   PATH.
#'
#' @param threads number of CPU threads for BWA and samtools.
#'
#' @param mem amount of RAM in GB (currently reserved; the streaming filter
#'   does not sort, so it needs no sort buffer).
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated before processing. FALSE resumes and skips only the lanes that
#'   already have complete output.
#'
#' @param overwrite.reference logical; if TRUE an existing BWA index in
#'   ref-index/ is deleted and rebuilt. The index is also rebuilt on its own
#'   when the contaminant reference files have changed. The index is removed
#'   when the function exits.
#'
#' @param quiet logical; if TRUE BWA and samtools stdout/stderr are suppressed.
#'
#' @return invisibly returns the summary data frame; writes decontaminated
#'   fastq.gz files to output.directory, per-lane contaminant counts to
#'   logs/sample_logs, and CSV summaries to logs/removeContamination_summary.csv
#'   and logs/removeContamination_contaminants.csv.
#'
#' @export

removeContamination = function(input.reads = "cleaned-reads",
                               output.directory = "decontaminated-reads",
                               decontamination.path = NULL,
                               map.match = 0.90,
                               samtools.path = NULL,
                               bwa.path = NULL,
                               threads = 1,
                               mem = 1,
                               overwrite = FALSE,
                               overwrite.reference = FALSE,
                               quiet = TRUE) {

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (file.exists(input.reads) == F){ stop("Input reads not found.") }
  if (is.null(decontamination.path) == TRUE){ stop("Please provide decontamination genomes / sequences.") }
  if (dir.exists(decontamination.path) == F){ stop("Decontamination directory not found.") }
  if (length(map.match) != 1 || is.numeric(map.match) == FALSE ||
      is.finite(map.match) == FALSE || map.match < 0 || map.match > 1){
    stop("map.match must be a number between 0 and 1.")
  }
  if (length(threads) != 1 || is.numeric(threads) == FALSE ||
      is.finite(threads) == FALSE || threads < 1) {
    stop("threads must be one positive number.")
  }
  if (length(overwrite) != 1 || is.logical(overwrite) == FALSE || is.na(overwrite) ||
      length(overwrite.reference) != 1 || is.logical(overwrite.reference) == FALSE ||
      is.na(overwrite.reference)) {
    stop("overwrite and overwrite.reference must be TRUE or FALSE.")
  }
  .checkDirectoryOverlap(input.reads, output.directory)

  #Checks that both programs are installed before any sample is processed
  bwa.command = .toolCommand("bwa", bwa.path)
  samtools.command = .toolCommand("samtools", samtools.path)

  # The combined reference index is temporary and is removed when this step exits.
  on.exit(unlink("ref-index", recursive = TRUE), add = TRUE)

  #Sets up the output directory
  if (dir.exists(output.directory) == F){
    dir.create(output.directory, recursive = TRUE)
  } else {
    if (overwrite == TRUE){ .resetDirectory(output.directory) }
  }#end else

  #Creates output directory
  if (dir.exists("logs/sample_logs") == F){ dir.create("logs/sample_logs", recursive = TRUE) }

  #Read in sample data
  input.reads = sub("/+$", "", input.reads)
  reads = .listFastqFiles(input.reads)
  read.names = .relativePaths(reads, input.reads)
  sample.names = .listSampleNames(input.reads)

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  #################################################
  ### Part A: build the combined contaminant index
  #################################################
  contig.map = .buildContaminantIndex(decontamination.path = decontamination.path,
                                      bwa.command = bwa.command,
                                      overwrite.reference = overwrite.reference,
                                      quiet = quiet)

  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             Lane = as.character(),
                             Task = as.character(),
                             Program = as.character(),
                             startPairs = as.numeric(),
                             removePairs = as.numeric(),
                             endPairs = as.numeric())

  contam.summary = data.frame(Sample = as.character(),
                              Lane = as.character(),
                              Contaminant = as.character(),
                              Accession = as.character(),
                              Reads = as.numeric())
  completed.samples = character()

  #Converts the identity threshold to the maximum allowed mismatch rate
  max.mismatch = 1 - map.match

  #Writes the awk program that splits contaminant pairs from clean pairs
  awk.script = tempfile(fileext = ".awk")
  on.exit(unlink(awk.script), add = TRUE)
  .writeContaminantAwk(awk.script)

  #Runs through each sample
  for (i in seq_along(sample.names)) {
    #################################################
    ### Part B: prepare for loading and checks
    #################################################
    sample.reads = .matchPrefix(reads, read.names, sample.names[i])

    #Returns a warning if reads are not found
    if (length(sample.reads) == 0 ){
      warning(sample.names[i], " does not have any reads present. Skipping.")
      next
    } #end if statement

    #Creates new directory
    out.path = paste0(output.directory, "/", sample.names[i])
    report.path = paste0("logs/sample_logs/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path, recursive = TRUE) }
    if (file.exists(report.path) == FALSE) { dir.create(report.path, recursive = TRUE) }

    lane.prefixes = .stripReadSuffix(sample.reads)

    for (j in seq_along(lane.prefixes)){

      lane.reads = .matchPrefix(reads, reads, lane.prefixes[j])

      #Returns a warning if reads are not found
      if (length(lane.reads) != 2 ){
        warning(lane.prefixes[j], " does not have a read pair present. Skipping.")
        next
      } #end if statement
      lane.reads = tryCatch(.orderReadPair(lane.reads), error = function(e) {
        warning(conditionMessage(e))
        return(NULL)
      })
      if (is.null(lane.reads) == TRUE) { next }

      lane.name = basename(lane.prefixes[j])
      lane.tag = gsub(".*_", "", lane.name)

      #sets up output reads
      outreads = c(paste0(out.path, "/", lane.name, "_READ1.fastq.gz"),
                   paste0(out.path, "/", lane.name, "_READ2.fastq.gz"))
      contam.csv = paste0(report.path, "/", lane.name, "_contamination-read-counts.csv")
      lane.csv = paste0(report.path, "/", lane.name, "_decontamination-summary.csv")
      counts.file = paste0(out.path, "/", lane.name, "_decontam-counts.txt")
      stats.file = paste0(out.path, "/", lane.name, "_decontam-stats.txt")
      metadata.file = paste0(report.path, "/", lane.name, "_decontamination-metadata.csv")
      metadata = .laneMetadata(lane.reads,
                               list(map.match = map.match,
                                    reference = attr(contig.map, "reference.identity")))

      # Reuses a lane only when its outputs, reports, inputs, reference, and
      # threshold all match the saved completion metadata.
      if (overwrite == FALSE && .laneComplete(outreads, require.size = FALSE) == TRUE &&
          file.exists(lane.csv) == TRUE && file.exists(contam.csv) == TRUE &&
          .metadataMatches(metadata.file, metadata) == TRUE) {
        lane.summary = tryCatch(read.csv(lane.csv, stringsAsFactors = FALSE),
                                error = function(e) NULL)
        lane.contam = tryCatch(read.csv(contam.csv, stringsAsFactors = FALSE),
                               error = function(e) NULL)
        if (is.null(lane.summary) == FALSE && is.null(lane.contam) == FALSE &&
            identical(names(lane.summary), names(summary.data)) &&
            identical(names(lane.contam), names(contam.summary))) {
          summary.data = rbind(summary.data, lane.summary)
          contam.summary = rbind(contam.summary, lane.contam)
          completed.samples = unique(c(completed.samples, sample.names[i]))
          print(paste0(lane.name, " is already complete. Skipping."))
          next
        }
      }
      if (overwrite == FALSE && .metadataConflicts(metadata.file, metadata) == TRUE) {
        stop(lane.name, " was completed with different inputs, reference, or threshold. ",
             "Use overwrite = TRUE to replace it.")
      }

      temp.outreads = vapply(outreads, function(output.file) {
        tempfile(pattern = paste0(basename(output.file), "-"),
                 tmpdir = out.path, fileext = ".fastq.gz")
      }, character(1))
      counts.file = tempfile(pattern = paste0(lane.name, "-counts-"), tmpdir = out.path)
      stats.file = tempfile(pattern = paste0(lane.name, "-stats-"), tmpdir = out.path)
      temp.lane.csv = tempfile(pattern = paste0(lane.name, "-summary-"),
                               tmpdir = report.path, fileext = ".csv")
      temp.contam.csv = tempfile(pattern = paste0(lane.name, "-contaminants-"),
                                 tmpdir = report.path, fileext = ".csv")
      on.exit(unlink(c(temp.outreads, counts.file, stats.file,
                       temp.lane.csv, temp.contam.csv)), add = TRUE)

      #################################################
      ### Part C: map, filter, and write the clean reads
      #################################################
      # One streaming pass does all of it. bwa writes the alignments, samtools
      # collate puts both mates of a pair together, the awk program splits the
      # contaminant pairs from the clean pairs, and samtools fastq writes the
      # clean pairs. Nothing is sorted by coordinate and no BAM is kept.
      .runPipeline(paste0(bwa.command, " mem -M -t ", threads, " ref-index/reference ",
                          shQuote(lane.reads[1]), " ", shQuote(lane.reads[2]),
                          " | ", samtools.command, " collate -O -u -@ ", threads, " - ",
                          " | ", samtools.command, " view -h -F 0x900 - ",
                          " | awk -v MAXMM=", max.mismatch,
                          " -v COUNTS=", shQuote(counts.file),
                          " -v STATS=", shQuote(stats.file),
                          " -f ", shQuote(awk.script), " - ",
                          " | ", samtools.command, " fastq -@ ", threads, " -n",
                          " -1 ", shQuote(temp.outreads[1]), " -2 ", shQuote(temp.outreads[2]),
                          " -0 /dev/null -s /dev/null -"),
                   quiet = quiet, task = "contaminant removal")

      #################################################
      ### Part D: gather the statistics
      #################################################
      # The awk program reports the pair totals, so no pass over the fastq files
      # is needed to count the reads.
      lane.stats = read.table(stats.file, sep = "\t", header = FALSE,
                              col.names = c("removed", "kept"))
      removed.pairs = lane.stats$removed[1]
      end.pairs = lane.stats$kept[1]

      contam.data = .readContaminantCounts(counts.file, contig.map)

      temp.remove = data.frame(Sample = sample.names[i],
                               Lane = lane.tag,
                               Task = "decontamination",
                               Program = "bwa",
                               startPairs = removed.pairs + end.pairs,
                               removePairs = removed.pairs,
                               endPairs = end.pairs)

      summary.data = rbind(summary.data, temp.remove)
      write.csv(temp.remove, file = temp.lane.csv, row.names = FALSE)

      if (nrow(contam.data) == 0) {
        lane.contam = contam.summary[0, , drop = FALSE]
      } else {
        lane.contam = data.frame(Sample = rep(sample.names[i], nrow(contam.data)),
                                 Lane = rep(lane.tag, nrow(contam.data)),
                                 Contaminant = contam.data$Contaminant,
                                 Accession = contam.data$Accession,
                                 Reads = contam.data$Reads,
                                 stringsAsFactors = FALSE)
      }

      contam.summary = rbind(contam.summary, lane.contam)
      write.csv(lane.contam, file = temp.contam.csv, row.names = FALSE)

      .publishFiles(c(temp.outreads, temp.lane.csv, temp.contam.csv),
                    c(outreads, lane.csv, contam.csv))
      .writeLaneMetadata(metadata, metadata.file)
      completed.samples = unique(c(completed.samples, sample.names[i]))

      unlink(c(counts.file, stats.file))

      print(paste0(lane.name, " completed decontamination!"))

    }#end lane j loop

    print(paste0(sample.names[i], " Completed decontamination removal!"))

  }#end sample i loop

  .appendSummary(summary.data, "logs/removeContamination_summary.csv",
                 replace.samples = completed.samples)
  .appendSummary(contam.summary, "logs/removeContamination_contaminants.csv",
                 replace.samples = completed.samples)

  return(invisible(summary.data))
}#end function


# Internal helper: builds the combined contaminant reference and its BWA index
# in ref-index/, and returns the contig to genome table. The index is rebuilt
# when the contaminant files have changed, so an old index cannot be used with a
# new contaminant list.
.buildContaminantIndex = function(decontamination.path = NULL,
                                  bwa.command = NULL,
                                  overwrite.reference = FALSE,
                                  quiet = TRUE) {

  active.file = file.path(decontamination.path, "active-references.csv")
  if (file.exists(active.file) == TRUE) {
    active.data = read.csv(active.file, stringsAsFactors = FALSE)
    if (!identical(names(active.data), "File")) {
      stop("The contaminant active-references.csv file must contain one File column.")
    }
    if (any(basename(active.data$File) != active.data$File) || any(duplicated(active.data$File))) {
      stop("The contaminant active reference list contains unsafe or duplicate file names.")
    }
    reference.list = file.path(decontamination.path, active.data$File)
  } else {
    reference.list = sort(list.files(decontamination.path, full.names = TRUE))
    reference.list = reference.list[grep("\\.fna\\.gz$|\\.fa\\.gz$|\\.fasta\\.gz$|\\.fna$|\\.fa$|\\.fasta$", reference.list)]
  }

  if (length(reference.list) == 0){
    stop("No contaminant reference files were found in ", decontamination.path, ".")
  }
  if (all(file.exists(reference.list)) == FALSE) {
    stop("The active contaminant reference list includes missing files: ",
         paste(basename(reference.list[!file.exists(reference.list)]), collapse = ", "))
  }
  for (reference.file in reference.list) { .checkFastaFile(reference.file) }

  # Checksums catch content changes even when a replacement has the same size
  # or an older modification time.
  manifest = paste0(normalizePath(reference.list), "\t",
                    unname(tools::md5sum(reference.list)))
  manifest.file = "ref-index/reference_files.txt"
  index.files = c("ref-index/reference.fa", "ref-index/reference.amb",
                  "ref-index/reference.ann", "ref-index/reference.bwt",
                  "ref-index/reference.pac", "ref-index/reference.sa",
                  "ref-index/contig_mapping.csv", manifest.file)

  index.stale = all(file.exists(index.files)) == FALSE
  if (index.stale == FALSE) {
    index.stale = !identical(readLines(manifest.file, warn = FALSE), manifest)
  }
  if (index.stale == TRUE && dir.exists("ref-index") == TRUE) {
    print("The contaminant reference files or index changed. Rebuilding the BWA index.")
  }

  if (overwrite.reference == TRUE || index.stale == TRUE) {
    temp.index = tempfile(pattern = "ref-index-", tmpdir = ".")
    dir.create(temp.index)
    on.exit(unlink(temp.index, recursive = TRUE), add = TRUE)
    index.prefix = file.path(temp.index, "reference")

    #Build contig-to-genome mapping before concatenating
    contig.map = data.frame(Contig = character(), Genome = character(), Accession = character())
    for (ref.file in reference.list) {
      file.base = gsub("\\.gz$", "", basename(ref.file))
      file.base = gsub("\\.fna$|\\.fa$|\\.fasta$", "", file.base)
      if (grepl("-", file.base)) {
        parts = strsplit(file.base, "-")[[1]]
        genome.name = parts[1]
        accession = paste(parts[-1], collapse = "-")
      } else {
        genome.name = file.base
        accession = file.base
      }
      headers = .runCommandOutput(paste0("gzip -cdf ", shQuote(ref.file), " | grep '^>'"),
                                  task = "contaminant header reading")
      contig.names = gsub("^>([^ ]+).*", "\\1", headers)
      contig.map = rbind(contig.map, data.frame(Contig = contig.names, Genome = genome.name, Accession = accession))
    }
    write.csv(contig.map, file = file.path(temp.index, "contig_mapping.csv"), row.names = FALSE)

    .runCommand(paste0("gzip -cdf ", paste(shQuote(reference.list), collapse = " "),
                       " > ", shQuote(paste0(index.prefix, ".fa"))),
                quiet = quiet, task = "contaminant reference concatenation", keep.stdout = TRUE)

    .runCommand(paste0(bwa.command, " index -p ", shQuote(index.prefix), " ",
                       shQuote(paste0(index.prefix, ".fa"))),
                quiet = quiet, task = "bwa index")

    temp.index.files = file.path(temp.index,
                                 c("reference.fa", "reference.amb", "reference.ann",
                                   "reference.bwt", "reference.pac", "reference.sa",
                                   "contig_mapping.csv"))
    if (all(file.exists(temp.index.files)) == FALSE) {
      stop("BWA did not create a complete contaminant index.")
    }
    writeLines(manifest, file.path(temp.index, "reference_files.txt"))

    old.index = NULL
    if (dir.exists("ref-index") == TRUE) {
      old.index = tempfile(pattern = "ref-index-old-", tmpdir = ".")
      if (file.rename("ref-index", old.index) == FALSE) {
        stop("Could not move the old contaminant index before replacement.")
      }
    }
    if (file.rename(temp.index, "ref-index") == FALSE) {
      if (is.null(old.index) == FALSE) { file.rename(old.index, "ref-index") }
      stop("Could not publish the completed contaminant index.")
    }
    if (is.null(old.index) == FALSE) { unlink(old.index, recursive = TRUE) }
  }#end rebuild index if

  contig.map = read.csv("ref-index/contig_mapping.csv", stringsAsFactors = FALSE)
  attr(contig.map, "reference.identity") = paste(manifest, collapse = ";")
  return(contig.map)
}#end .buildContaminantIndex


# Internal helper: turns the per-contig alignment counts into per-genome counts.
.readContaminantCounts = function(counts.file = NULL,
                                  contig.map = NULL) {

  empty.result = data.frame(Contaminant = character(), Accession = character(),
                            Reads = numeric(), stringsAsFactors = FALSE)

  if (file.exists(counts.file) == FALSE) { return(empty.result) }
  if (file.info(counts.file)$size == 0) { return(empty.result) }

  table.dat = read.table(counts.file, sep = "\t", header = FALSE,
                         col.names = c("Organism", "Count"), stringsAsFactors = FALSE)
  table.dat = table.dat[table.dat$Organism != "*", , drop = FALSE]
  if (nrow(table.dat) == 0) { return(empty.result) }

  table.dat = merge(table.dat, contig.map, by.x = "Organism", by.y = "Contig", all.x = TRUE)
  table.dat$Genome[is.na(table.dat$Genome)] = table.dat$Organism[is.na(table.dat$Genome)]
  table.dat$Accession[is.na(table.dat$Accession)] = table.dat$Organism[is.na(table.dat$Accession)]

  contam.data = aggregate(table.dat$Count, FUN = sum,
                          by = list(table.dat$Genome, table.dat$Accession))
  colnames(contam.data) = c("Contaminant", "Accession", "Reads")

  return(contam.data)
}#end .readContaminantCounts
