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
#'   row and falls back to the ENA file layout when the column is absent. A
#'   requested filter requires the LibraryLayout column. The ENA file report
#'   still verifies the resolved file layout.
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

sraDownload = function(sra.info.file = NULL,
                       sample.name.column = NULL,
                       output.directory = NULL,
                       filter.library.strategy = NULL,
                       filter.library.layout = NULL,
                       max.retries = 3,
                       retry.delay = 10,
                       skip.not.found = TRUE,
                       overwrite = FALSE,
                       quiet = FALSE) {

  #Quick checks
  if (is.null(sra.info.file) == TRUE){ stop("Please provide an sra.info.file path.") }
  if (file.exists(sra.info.file) == FALSE){ stop("sra.info.file not found: ", sra.info.file) }
  if (is.null(output.directory) == TRUE){ stop("Please provide an output.directory.") }
  if (length(max.retries) != 1 || is.numeric(max.retries) == FALSE ||
      is.finite(max.retries) == FALSE || max.retries < 1 || max.retries %% 1 != 0) {
    stop("max.retries must be one positive integer.")
  }
  if (length(retry.delay) != 1 || is.numeric(retry.delay) == FALSE ||
      is.finite(retry.delay) == FALSE || retry.delay < 0) {
    stop("retry.delay must be one non-negative number.")
  }
  if (length(skip.not.found) != 1 || is.logical(skip.not.found) == FALSE || is.na(skip.not.found)) {
    stop("skip.not.found must be TRUE or FALSE.")
  }
  if (length(overwrite) != 1 || is.logical(overwrite) == FALSE || is.na(overwrite)) {
    stop("overwrite must be TRUE or FALSE.")
  }

  #Reads the SRA information table
  sra.data = read.csv(sra.info.file, stringsAsFactors = FALSE)

  if ("Run" %in% names(sra.data) == FALSE){
    stop("sra.info.file must contain a 'Run' column with SRR/ERR/DRR accessions.")
  }
  sra.data$Run = trimws(as.character(sra.data$Run))
  if (any(is.na(sra.data$Run)) || any(nchar(sra.data$Run) == 0)) {
    stop("The Run column cannot contain missing or blank accessions.")
  }

  # Optional row filters
  if (is.null(filter.library.strategy) == FALSE &&
      "LibraryStrategy" %in% names(sra.data) == FALSE){
    stop("filter.library.strategy was requested, but LibraryStrategy is absent from sra.info.file.")
  }
  if (is.null(filter.library.strategy) == FALSE){
    sra.data = sra.data[which(sra.data$LibraryStrategy == filter.library.strategy), , drop = FALSE]
  }

  if (is.null(filter.library.layout) == FALSE &&
      "LibraryLayout" %in% names(sra.data) == FALSE){
    stop("filter.library.layout was requested, but LibraryLayout is absent from sra.info.file.")
  }
  if (is.null(filter.library.layout) == FALSE){
    filter.library.layout = toupper(filter.library.layout)
    sra.data = sra.data[which(toupper(sra.data$LibraryLayout) == filter.library.layout), , drop = FALSE]
  }

  if (nrow(sra.data) == 0){
    stop("No rows remain in sra.info.file after applying filters.")
  }

  .checkFileOutsideOutput(sra.info.file, output.directory)

  #Sets up the output directory
  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) { .resetDirectory(output.directory) }
  } else {
    dir.create(output.directory, recursive = TRUE)
  }
  if (dir.exists("logs/sample_logs") == FALSE) { dir.create("logs/sample_logs", recursive = TRUE) }

  #Builds the sample names
  # Priority:
  #   1. sample.name.column explicitly set -> use that column directly
  #   2. ScientificName + SampleName both present -> Genus_species_SampleName
  #      (falls back to Genus_species_Run for rows where SampleName is blank)
  #   3. ScientificName only              -> Genus_species_Run
  #   4. Neither                          -> Run accession alone
  if (is.null(sample.name.column) == FALSE) {
    if (sample.name.column %in% names(sra.data) == FALSE){
      stop("sample.name.column '", sample.name.column, "' not found in sra.info.file.")
    }
    sra.data$sample.name = as.character(sra.data[[sample.name.column]])
    if (any(is.na(sra.data$sample.name)) ||
        any(nchar(trimws(sra.data$sample.name)) == 0)) {
      stop("sample.name.column cannot contain missing or blank values.")
    }
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
  original.names = trimws(as.character(sra.data$sample.name))
  blank.names = is.na(original.names) | nchar(original.names) == 0
  original.names[blank.names] = sra.data$Run[blank.names]
  sra.data$sample.name = .sanitizeName(original.names)
  name.map = unique(data.frame(original = original.names,
                               clean = sra.data$sample.name,
                               stringsAsFactors = FALSE))
  if (any(duplicated(name.map$clean))) {
    stop("Distinct SRA sample names become identical after filename cleaning.")
  }

  #Finds the files for one accession
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

  #Runs the download
  # Multiple SRR accessions that resolve to the same sample name (same
  # ScientificName + SampleName) are treated as sequencing lanes of a single
  # individual, exactly as dropboxDownload handles multi-lane samples.
  # They are downloaded as L001, L002, ... and share one Sample entry in the
  # rename CSV so organizeReads merges them automatically.
  unique.samples = unique(sra.data$sample.name)
  n.total = length(unique.samples)
  rename.out = data.frame(File = character(), Sample = character(),
                          stringsAsFactors = FALSE)
  excluded.runs = data.frame(Run = character(), Sample = character(),
                             Reason = character(), stringsAsFactors = FALSE)

  for (i in seq_len(n.total)) {

    samp = unique.samples[i]
    samp.rows = sra.data[sra.data$sample.name == samp, ]
    n.lanes = nrow(samp.rows)

    if (quiet == FALSE) { message(sprintf("[%d/%d] %s  (%d run(s))", i, n.total, samp, n.lanes)) }

    sample.log = file.path("logs/sample_logs", samp)
    dir.create(sample.log, recursive = TRUE, showWarnings = FALSE)
    sentinel = file.path(sample.log, paste0(samp, "_sra-metadata.csv"))
    sample.metadata = c("sample" = samp,
                        "accessions" = paste(samp.rows$Run, collapse = ";"),
                        "layouts" = if ("LibraryLayout" %in% names(samp.rows)) {
                          paste(samp.rows$LibraryLayout, collapse = ";")
                        } else {
                          "not supplied"
                        })
    if (overwrite == FALSE && .metadataConflicts(sentinel, sample.metadata) == TRUE) {
      stop(samp, " has a completed SRA download with a different accession list or order. ",
           "Use overwrite = TRUE to replace it.")
    }

    #Runs each lane
    all.lanes.ok = TRUE

    for (j in seq_len(n.lanes)) {

      acc = samp.rows$Run[j]
      if ("LibraryLayout" %in% names(samp.rows)) {
        layout = samp.rows$LibraryLayout[j]
        if (is.na(layout) || nchar(trimws(layout)) == 0) { layout = "PAIRED" }
      } else {
        layout = "PAIRED"
      }
      lane.tag = sprintf("L%03d", j)

      if (quiet == FALSE && n.lanes > 1){
        message(sprintf("  lane %d/%d (%s)", j, n.lanes, acc))
      }

      run.files = .runFiles(acc, layout)
      is.paired = is.na(run.files$r2) == FALSE

      if (identical(filter.library.layout, "PAIRED") && is.paired == FALSE) {
        msg = paste0(acc, " resolved to one FASTQ file, but paired reads were requested.")
        excluded.runs = rbind(excluded.runs,
                              data.frame(Run = acc, Sample = samp, Reason = msg,
                                         stringsAsFactors = FALSE))
        if (skip.not.found == FALSE) { stop(msg) }
        warning(msg)
        all.lanes.ok = FALSE
        next
      }

      # Destination paths for this lane
      r1.dest = file.path(output.directory, paste0(samp, "_", lane.tag, "_READ1.fastq.gz"))
      r2.dest = if (is.paired)
                  file.path(output.directory, paste0(samp, "_", lane.tag, "_READ2.fastq.gz"))
                else NULL

      metadata.file = file.path(sample.log,
                                paste0(samp, "_", lane.tag, "_sra-lane-metadata.csv"))
      metadata = c("accession" = acc, "layout" = if (is.paired) "PAIRED" else "SINGLE",
                   "source.1" = run.files$r1,
                   "source.2" = if (is.paired) run.files$r2 else "",
                   "md5.1" = run.files$md5.r1,
                   "md5.2" = if (is.paired) run.files$md5.r2 else "")

      # Skip this lane only when its files and matching completion metadata exist.
      if (.laneComplete(c(r1.dest, r2.dest), metadata.file = metadata.file,
                        metadata = metadata) == TRUE) {
        if (quiet == FALSE) { message("    ", lane.tag, " files exist -- skipping") }
        rename.out = rbind(rename.out,
                           data.frame(File = paste0(samp, "_", lane.tag),
                                      Sample = samp, stringsAsFactors = FALSE))
        next
      }
      if (overwrite == FALSE && .metadataConflicts(metadata.file, metadata) == TRUE) {
        stop(samp, " ", lane.tag, " was downloaded from a different SRA run. ",
             "Use overwrite = TRUE to replace it.")
      }

      output.files = c(r1.dest, r2.dest)
      temp.files = vapply(output.files, function(output.file) {
        tempfile(pattern = paste0(basename(output.file), "-"),
                 tmpdir = output.directory, fileext = ".fastq.gz")
      }, character(1))
      on.exit(unlink(temp.files), add = TRUE)

      # Download READ1
      r1.ok = .dl(run.files$r1, temp.files[1], max.retries, retry.delay, quiet,
                  run.files$md5.r1)
      if (!r1.ok) {
        msg = sprintf("  READ1 download failed for %s (%s) after %d attempts",
                      acc, lane.tag, max.retries)
        excluded.runs = rbind(excluded.runs,
                              data.frame(Run = acc, Sample = samp, Reason = msg,
                                         stringsAsFactors = FALSE))
        if (skip.not.found == FALSE) {
          stop(msg)
        } else {
          warning(msg)
          all.lanes.ok = FALSE
          next
        }
      }

      # Download READ2 (PAIRED only)
      if (is.paired) {
        r2.ok = .dl(run.files$r2, temp.files[2], max.retries, retry.delay, quiet,
                    run.files$md5.r2)
        if (!r2.ok) {
          unlink(temp.files)
          msg = sprintf("  READ2 download failed for %s (%s) after %d attempts",
                        acc, lane.tag, max.retries)
          excluded.runs = rbind(excluded.runs,
                                data.frame(Run = acc, Sample = samp, Reason = msg,
                                           stringsAsFactors = FALSE))
          if (skip.not.found == FALSE) {
            stop(msg)
          } else {
            warning(msg)
            all.lanes.ok = FALSE
            next
          }
        }
      }

      .publishFiles(temp.files, output.files)
      .writeLaneMetadata(metadata, metadata.file)

      if (quiet == FALSE) { message("    ", lane.tag, " done") }
      rename.out = rbind(rename.out,
                         data.frame(File = paste0(samp, "_", lane.tag),
                                    Sample = samp, stringsAsFactors = FALSE))
    } # end lane loop

    # Write sample metadata only when every requested lane succeeded.
    if (all.lanes.ok) {
      .writeLaneMetadata(sample.metadata, sentinel)
    }

  } # end sample loop

  #Writes the rename CSV
  # The CSV is quoted so that a sample name with a comma cannot break the table
  write.csv(rename.out,
            file = "file_rename_sra.csv",
            row.names = FALSE)
  if (nrow(excluded.runs) > 0) {
    dir.create("logs", showWarnings = FALSE)
    write.csv(excluded.runs, "logs/sraDownload_excluded-runs.csv", row.names = FALSE)
  }

  if (quiet == FALSE){
    message("\nDone. ", nrow(rename.out), " sample(s) recorded in file_rename_sra.csv.")
  }

  return(invisible(rename.out))

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
