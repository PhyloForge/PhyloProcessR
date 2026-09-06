#' @title sraDownload
#'
#' @description Downloads paired-end or single-end FASTQ files from NCBI's
#'   Sequence Read Archive (SRA) using the European Nucleotide Archive (ENA)
#'   HTTPS mirrors, which provide pre-formatted FASTQ.gz files without
#'   requiring any external SRA toolkit installation. Accepts an SraRunInfo
#'   CSV file exported from the NCBI SRA Run Selector, or any CSV that
#'   contains at minimum a 'Run' column with SRR/ERR/DRR accession numbers.
#'   The ENA file report gives the exact file paths and their MD5 checksums,
#'   so runs with an unusual file layout are handled and every download is
#'   verified. Downloaded files are named to the standard PhyloProcessR
#'   convention (SampleName_L001_READ1/2.fastq.gz) and a file_rename_sra.csv is
#'   written in the working directory for direct use with organizeReads.
#'
#' @param sra.info.file character; path to the SraRunInfo CSV. Must contain
#'   at minimum a 'Run' column. The full NCBI SRA Run Selector export (all
#'   columns) is accepted directly; only the columns described below are used.
#'
#' @param sample.name.column character or NULL; name of a column in
#'   sra.info.file to use directly as sample names. If NULL (default), names
#'   are built automatically in priority order: (1) if both 'ScientificName'
#'   and 'SampleName' columns are present the name is Genus_species_SampleName
#'   (e.g. Hylarana_macrodactyla_CAS12345), falling back to
#'   Genus_species_SRRaccession for any row where SampleName is blank; (2) if
#'   only 'ScientificName' is present, Genus_species_SRRaccession; (3)
#'   otherwise the bare SRR accession alone. Every name is cleaned so that only
#'   letters, digits, and the characters . _ - remain.
#'
#' @param output.directory character; local directory where the FASTQ.gz files
#'   will be saved. Created if it does not exist.
#'
#' @param filter.library.strategy character or NULL; if provided, only rows
#'   whose 'LibraryStrategy' column matches this string are downloaded (e.g.
#'   "Targeted-Capture"). NULL (default) downloads all rows.
#'
#' @param filter.library.layout character or NULL; restrict downloads to
#'   "PAIRED" or "SINGLE". NULL (default) reads the LibraryLayout column per
#'   row; falls back to PAIRED if the column is absent. The ENA file report
#'   overrides this when it is available.
#'
#' @param max.retries integer; number of download attempts per file before
#'   giving up. Default 3.
#'
#' @param retry.delay numeric; seconds to wait between retry attempts.
#'   Default 10.
#'
#' @param skip.not.found logical; if TRUE a warning is printed and the sample
#'   is skipped when all retry attempts fail. If FALSE an error is raised.
#'   Default TRUE.
#'
#' @param overwrite logical; if TRUE the output directory is deleted and
#'   recreated and every sample is re-downloaded from scratch. Default FALSE.
#'
#' @param quiet logical; suppress progress messages. Default FALSE.
#'
#' @return invisibly returns a data.frame with File and Sample columns (the
#'   same content written to file_rename_sra.csv). Side effects: FASTQ.gz
#'   files in output.directory and file_rename_sra.csv in the working
#'   directory.
#'
#' @export

sraDownload = function(sra.info.file           = NULL,
                       sample.name.column       = NULL,
                       output.directory         = NULL,
                       filter.library.strategy  = NULL,
                       filter.library.layout    = NULL,
                       max.retries              = 3,
                       retry.delay              = 10,
                       skip.not.found           = TRUE,
                       overwrite                = FALSE,
                       quiet                    = FALSE) {

  # -- Argument checks ---------------------------------------------------------
  if (is.null(sra.info.file))      stop("Please provide an sra.info.file path.")
  if (!file.exists(sra.info.file)) stop("sra.info.file not found: ", sra.info.file)
  if (is.null(output.directory))   stop("Please provide an output.directory.")

  # -- Output directory ---------------------------------------------------------
  if (dir.exists(output.directory)) {
    if (overwrite) { .resetDirectory(output.directory) }
  } else {
    dir.create(output.directory, recursive = TRUE)
  }

  # -- Read SRA info table ------------------------------------------------------
  sra.data = read.csv(sra.info.file, stringsAsFactors = FALSE)

  if (!"Run" %in% names(sra.data))
    stop("sra.info.file must contain a 'Run' column with SRR/ERR/DRR accessions.")

  # Optional row filters
  if (!is.null(filter.library.strategy) && "LibraryStrategy" %in% names(sra.data))
    sra.data = sra.data[sra.data$LibraryStrategy == filter.library.strategy, ]

  if (!is.null(filter.library.layout) && "LibraryLayout" %in% names(sra.data))
    sra.data = sra.data[sra.data$LibraryLayout == filter.library.layout, ]

  if (nrow(sra.data) == 0)
    stop("No rows remain in sra.info.file after applying filters.")

  # -- Build sample names -------------------------------------------------------
  # Priority:
  #   1. sample.name.column explicitly set -> use that column directly
  #   2. ScientificName + SampleName both present -> Genus_species_SampleName
  #      (falls back to Genus_species_Run for rows where SampleName is blank)
  #   3. ScientificName only              -> Genus_species_Run
  #   4. Neither                          -> Run accession alone
  if (!is.null(sample.name.column)) {
    if (!sample.name.column %in% names(sra.data))
      stop("sample.name.column '", sample.name.column, "' not found in sra.info.file.")
    sra.data$sample.name = as.character(sra.data[[sample.name.column]])
  } else if ("ScientificName" %in% names(sra.data)) {
    sci = trimws(sra.data$ScientificName)
    if ("SampleName" %in% names(sra.data)) {
      sn = trimws(as.character(sra.data$SampleName))
      # Use SampleName where non-empty; fall back to Run accession for blank rows
      specimen.id = ifelse(nchar(sn) > 0 & !is.na(sn), sn, sra.data$Run)
      sra.data$sample.name = paste0(sci, "_", specimen.id)
    } else {
      sra.data$sample.name = paste0(sci, "_", sra.data$Run)
    }
  } else {
    sra.data$sample.name = sra.data$Run
  }

  # A sample name reaches a file path and a CSV field, so a space, a comma, or a
  # regular expression character is replaced here.
  sra.data$sample.name = .sanitizeName(sra.data$sample.name)
  blank.names = nchar(sra.data$sample.name) == 0
  if (any(blank.names)) { sra.data$sample.name[blank.names] = sra.data$Run[blank.names] }

  # -- Internal: file locations for an accession -------------------------------
  # The ENA file report gives the true file paths and their MD5 checksums. It
  # handles runs that hold only one file and runs that hold an extra unpaired
  # file. The URL pattern below is the fallback when the report is unavailable.
  .runFiles = function(acc, layout = "PAIRED") {

    report = .enaFileReport(acc)
    if (is.null(report) == FALSE) { return(report) }

    urls = .enaUrls(acc, layout)
    return(list(r1 = unname(urls["r1"]),
                r2 = if (length(urls) > 1) unname(urls["r2"]) else NA_character_,
                md5.r1 = NA_character_,
                md5.r2 = NA_character_))
  }

  # -- Main download loop -------------------------------------------------------
  # Multiple SRR accessions that resolve to the same sample name (same
  # ScientificName + SampleName) are treated as sequencing lanes of a single
  # individual, exactly as dropboxDownload handles multi-lane samples.
  # They are downloaded as L001, L002, ... and share one Sample entry in the
  # rename CSV so organizeReads merges them automatically.
  #
  # Sentinel: one file per sample named SampleName.fastq.sra_done, written
  # only after ALL lanes for that sample complete. The .fastq.* component
  # matches the gsub strip used by fastqStats / readStats so the sentinel
  # collapses to SampleName and is never treated as a separate sample.
  # The dot separator (not underscore) means the prefix match in those same
  # functions never picks it up as a read file.
  unique.samples  = unique(sra.data$sample.name)
  n.total         = length(unique.samples)
  rename.out      = data.frame(File = character(), Sample = character(),
                               stringsAsFactors = FALSE)

  for (i in seq_len(n.total)) {

    samp      = unique.samples[i]
    samp.rows = sra.data[sra.data$sample.name == samp, ]
    n.lanes   = nrow(samp.rows)

    if (!quiet) message(sprintf("[%d/%d] %s  (%d run(s))", i, n.total, samp, n.lanes))

    # Sample-level sentinel -- fast skip when all lanes completed in a prior run.
    sentinel = file.path(output.directory, paste0(samp, ".fastq.sra_done"))
    if (file.exists(sentinel)) {
      if (!quiet) message("  all lanes already completed -- skipping")
      # Recover rename entries from the files that actually exist on disk
      existing.lanes = list.files(output.directory)
      existing.lanes = existing.lanes[startsWith(existing.lanes, paste0(samp, "_")) &
                                        endsWith(existing.lanes, "_READ1.fastq.gz")]
      lane.tags = sub("_READ1\\.fastq\\.gz$", "", substring(existing.lanes, nchar(samp) + 2))
      for (lt in sort(lane.tags)) {
        rename.out = rbind(rename.out,
                           data.frame(File   = paste0(samp, "_", lt),
                                      Sample = samp, stringsAsFactors = FALSE))
      }
      next
    }

    # -- Inner lane loop --------------------------------------------------------
    all.lanes.ok = TRUE

    for (j in seq_len(n.lanes)) {

      acc       = samp.rows$Run[j]
      layout    = if ("LibraryLayout" %in% names(samp.rows)) samp.rows$LibraryLayout[j] else "PAIRED"
      lane.tag  = sprintf("L%03d", j)

      if (!quiet && n.lanes > 1)
        message(sprintf("  lane %d/%d (%s)", j, n.lanes, acc))

      run.files = .runFiles(acc, layout)
      is.paired = is.na(run.files$r2) == FALSE

      # Destination paths for this lane
      r1.dest = file.path(output.directory, paste0(samp, "_", lane.tag, "_READ1.fastq.gz"))
      r2.dest = if (is.paired)
                  file.path(output.directory, paste0(samp, "_", lane.tag, "_READ2.fastq.gz"))
                else NULL

      # Skip this lane if its files already exist (prior partial run)
      if (.laneComplete(c(r1.dest, r2.dest)) == TRUE) {
        if (!quiet) message("    ", lane.tag, " files exist -- skipping")
        rename.out = rbind(rename.out,
                           data.frame(File   = paste0(samp, "_", lane.tag),
                                      Sample = samp, stringsAsFactors = FALSE))
        next
      }

      # Download READ1
      r1.ok = .dl(run.files$r1, r1.dest, max.retries, retry.delay, quiet, run.files$md5.r1)
      if (!r1.ok) {
        if (file.exists(r1.dest)) file.remove(r1.dest)
        msg = sprintf("  READ1 download failed for %s (%s) after %d attempts",
                      acc, lane.tag, max.retries)
        if (!skip.not.found) stop(msg) else { warning(msg); all.lanes.ok = FALSE; next }
      }

      # Download READ2 (PAIRED only)
      if (is.paired) {
        r2.ok = .dl(run.files$r2, r2.dest, max.retries, retry.delay, quiet, run.files$md5.r2)
        if (!r2.ok) {
          if (file.exists(r1.dest)) file.remove(r1.dest)
          if (file.exists(r2.dest)) file.remove(r2.dest)
          msg = sprintf("  READ2 download failed for %s (%s) after %d attempts",
                        acc, lane.tag, max.retries)
          if (!skip.not.found) stop(msg) else { warning(msg); all.lanes.ok = FALSE; next }
        }
      }

      if (!quiet) message("    ", lane.tag, " done")
      rename.out = rbind(rename.out,
                         data.frame(File   = paste0(samp, "_", lane.tag),
                                    Sample = samp, stringsAsFactors = FALSE))
    } # end lane loop

    # Write sample-level sentinel only when every lane succeeded
    if (all.lanes.ok) {
      writeLines(c(
        paste0("sample:    ", samp),
        paste0("lanes:     ", n.lanes),
        paste0("accessions:", paste(samp.rows$Run, collapse = " ")),
        paste0("completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      ), sentinel)
    }

  } # end sample loop

  # -- Write rename CSV ---------------------------------------------------------
  # The CSV is quoted so that a sample name with a comma cannot break the table
  write.csv(rename.out,
            file      = "file_rename_sra.csv",
            row.names = FALSE)

  if (!quiet)
    message("\nDone. ", nrow(rename.out), " sample(s) recorded in file_rename_sra.csv.")

  invisible(rename.out)

} # end sraDownload


# Internal helper: asks the ENA file report for the FASTQ paths and MD5 sums of
# one run accession. Returns NULL when the report is unavailable, so the caller
# can fall back to the URL pattern.
.enaFileReport = function(acc = NULL) {

  report.url = paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", acc,
                      "&result=read_run&fields=fastq_ftp,fastq_md5&format=tsv")

  report = tryCatch(utils::read.delim(report.url, stringsAsFactors = FALSE),
                    error = function(e) NULL, warning = function(w) NULL)

  if (is.null(report) == TRUE) { return(NULL) }
  if (nrow(report) == 0) { return(NULL) }
  if ("fastq_ftp" %in% names(report) == FALSE) { return(NULL) }
  if (is.na(report$fastq_ftp[1]) || nchar(report$fastq_ftp[1]) == 0) { return(NULL) }

  file.paths = unlist(strsplit(report$fastq_ftp[1], ";", fixed = TRUE))
  file.md5 = rep(NA_character_, length(file.paths))
  if ("fastq_md5" %in% names(report) == TRUE) {
    md5.values = unlist(strsplit(report$fastq_md5[1], ";", fixed = TRUE))
    if (length(md5.values) == length(file.paths)) { file.md5 = md5.values }
  }

  file.names = basename(file.paths)
  file.urls = paste0("https://", sub("^https?://", "", file.paths))

  # A paired run holds _1 and _2 files. A single run holds one file, which ENA
  # may also add to a paired run to hold the orphan reads.
  first.mate = endsWith(file.names, "_1.fastq.gz")
  second.mate = endsWith(file.names, "_2.fastq.gz")

  if (any(first.mate) && any(second.mate)) {
    return(list(r1 = file.urls[first.mate][1],
                r2 = file.urls[second.mate][1],
                md5.r1 = file.md5[first.mate][1],
                md5.r2 = file.md5[second.mate][1]))
  }

  return(list(r1 = file.urls[1],
              r2 = NA_character_,
              md5.r1 = file.md5[1],
              md5.r2 = NA_character_))
}#end .enaFileReport


# Internal helper: builds the ENA HTTPS URLs for an accession from the standard
# path pattern. This is the fallback when the ENA file report is unavailable.
# URL structure:
#   ftp.sra.ebi.ac.uk/vol1/fastq/{first6}/[subdir]/{acc}/{acc}_[1|2].fastq.gz
# subdir is derived from the accession length:
#   <= 9 chars : no subdir
#   10 chars   : 00{last1}
#   11 chars   : 0{last2}
#   12 chars   : {last3}
.enaUrls = function(acc = NULL,
                    layout = "PAIRED") {

  n    = nchar(acc)
  f6   = substr(acc, 1, 6)
  sub.dir = if      (n <= 9)  ""
            else if (n == 10) paste0("/00", substr(acc, n,     n  ))
            else if (n == 11) paste0("/0",  substr(acc, n - 1, n  ))
            else               paste0("/",  substr(acc, n - 2, n  ))
  base = paste0("https://ftp.sra.ebi.ac.uk/vol1/fastq/", f6, sub.dir, "/", acc, "/")

  if (toupper(layout) == "PAIRED") {
    return(c(r1 = paste0(base, acc, "_1.fastq.gz"),
             r2 = paste0(base, acc, "_2.fastq.gz")))
  }

  return(c(r1 = paste0(base, acc, ".fastq.gz")))
}#end .enaUrls


# Internal helper: downloads one file with retries and verifies it.
# utils::download.file only *warns* on timeout or length mismatch -- it never
# throws an error -- so a plain tryCatch misses truncated files. We use
# withCallingHandlers to intercept those warnings and treat them as failures.
# The MD5 sum from the ENA file report gives a second check that catches a
# download that finished but holds the wrong bytes.
# The global timeout option is raised to 3600 s for the duration of the call
# (large FASTQ files easily exceed the 60-second default).
.dl = function(src.url, dest.path, max.retries, retry.delay, quiet, expected.md5 = NA_character_) {
  old.timeout = getOption("timeout")
  options(timeout = 3600)
  on.exit(options(timeout = old.timeout), add = TRUE)

  for (attempt in seq_len(max.retries)) {
    bad.warn = FALSE
    ok = tryCatch({
      withCallingHandlers(
        utils::download.file(src.url, dest.path, mode = "wb", quiet = TRUE),
        warning = function(w) {
          msg = conditionMessage(w)
          if (grepl("downloaded length|Timeout|timed out", msg, ignore.case = TRUE))
            bad.warn <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
      !bad.warn && file.exists(dest.path) && file.size(dest.path) > 0
    }, error = function(e) FALSE)

    if (ok && is.na(expected.md5) == FALSE && nchar(expected.md5) > 0) {
      file.md5 = unname(tools::md5sum(dest.path))
      if (identical(file.md5, expected.md5) == FALSE) {
        if (!quiet) message("    checksum did not match -- the file is incomplete")
        ok = FALSE
      }
    }

    if (ok) return(TRUE)

    if (file.exists(dest.path)) file.remove(dest.path)
    if (attempt < max.retries) {
      if (!quiet) message("    attempt ", attempt, " failed -- retrying in ", retry.delay, "s")
      Sys.sleep(retry.delay)
    }
  }
  FALSE
}#end .dl
