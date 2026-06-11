#' @title sraDownload
#'
#' @description Downloads paired-end or single-end FASTQ files from NCBI's
#'   Sequence Read Archive (SRA) using the European Nucleotide Archive (ENA)
#'   HTTPS mirrors, which provide pre-formatted FASTQ.gz files without
#'   requiring any external SRA toolkit installation. Accepts an SraRunInfo
#'   CSV file exported from the NCBI SRA Run Selector, or any CSV that
#'   contains at minimum a 'Run' column with SRR/ERR/DRR accession numbers.
#'   Downloaded files are named to the standard PhyloProcessR convention
#'   (SampleName_L001_READ1.fastq.gz / READ2.fastq.gz) and a
#'   file_rename_sra.csv is written in the working directory for direct use
#'   with organizeReads.
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
#'   otherwise the bare SRR accession alone.
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
#'   row; falls back to PAIRED if the column is absent.
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

  # ── Argument checks ─────────────────────────────────────────────────────────
  if (is.null(sra.info.file))      stop("Please provide an sra.info.file path.")
  if (!file.exists(sra.info.file)) stop("sra.info.file not found: ", sra.info.file)
  if (is.null(output.directory))   stop("Please provide an output.directory.")

  # ── Output directory ─────────────────────────────────────────────────────────
  if (dir.exists(output.directory)) {
    if (overwrite) {
      system(paste0("rm -r ", shQuote(output.directory)))
      dir.create(output.directory, recursive = TRUE)
    }
  } else {
    dir.create(output.directory, recursive = TRUE)
  }

  # ── Read SRA info table ──────────────────────────────────────────────────────
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

  # ── Build sample names ───────────────────────────────────────────────────────
  # Priority:
  #   1. sample.name.column explicitly set → use that column directly
  #   2. ScientificName + SampleName both present → Genus_species_SampleName
  #      (falls back to Genus_species_Run for rows where SampleName is blank)
  #   3. ScientificName only              → Genus_species_Run
  #   4. Neither                          → Run accession alone
  if (!is.null(sample.name.column)) {
    if (!sample.name.column %in% names(sra.data))
      stop("sample.name.column '", sample.name.column, "' not found in sra.info.file.")
    sra.data$sample.name = as.character(sra.data[[sample.name.column]])
  } else if ("ScientificName" %in% names(sra.data)) {
    sci = gsub("[[:space:]]+", "_", trimws(sra.data$ScientificName))
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

  # ── Internal: ENA HTTPS URL(s) for an accession ─────────────────────────────
  # ENA mirrors all public SRA data as pre-formatted FASTQ.gz.
  # URL structure:
  #   ftp.sra.ebi.ac.uk/vol1/fastq/{first6}/[subdir]/{acc}/{acc}_[1|2].fastq.gz
  # subdir is derived from the accession length:
  #   <= 9 chars : no subdir
  #   10 chars   : 00{last1}
  #   11 chars   : 0{last2}
  #   12 chars   : {last3}
  .ena.urls = function(acc, layout = "PAIRED") {
    n    = nchar(acc)
    f6   = substr(acc, 1, 6)
    sub  = if      (n <= 9)  ""
            else if (n == 10) paste0("/00", substr(acc, n,     n  ))
            else if (n == 11) paste0("/0",  substr(acc, n - 1, n  ))
            else               paste0("/",  substr(acc, n - 2, n  ))
    base = paste0("https://ftp.sra.ebi.ac.uk/vol1/fastq/", f6, sub, "/", acc, "/")
    if (toupper(layout) == "PAIRED") {
      c(r1 = paste0(base, acc, "_1.fastq.gz"),
        r2 = paste0(base, acc, "_2.fastq.gz"))
    } else {
      c(r1 = paste0(base, acc, ".fastq.gz"))
    }
  }

  # ── Internal: download one file with retries ─────────────────────────────────
  # utils::download.file only *warns* on timeout or length mismatch — it never
  # throws an error — so a plain tryCatch misses truncated files. We use
  # withCallingHandlers to intercept those warnings and treat them as failures.
  # The global timeout option is raised to 3600 s for the duration of the call
  # (large FASTQ files easily exceed the 60-second default).
  .dl = function(src.url, dest.path, max.retries, retry.delay, quiet) {
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

      if (ok) return(TRUE)

      if (file.exists(dest.path)) file.remove(dest.path)
      if (attempt < max.retries) {
        if (!quiet) message("    attempt ", attempt, " failed — retrying in ", retry.delay, "s")
        Sys.sleep(retry.delay)
      }
    }
    FALSE
  }

  # ── Main download loop ───────────────────────────────────────────────────────
  # Sentinels are stored in a hidden subdirectory so they are never picked up
  # by downstream tools (fastqStats, organizeReads, etc.) that scan the reads
  # directory for sample files.
  sentinels.dir = file.path(output.directory, ".sra-sentinels")
  dir.create(sentinels.dir, showWarnings = FALSE, recursive = TRUE)

  n.total    = nrow(sra.data)
  rename.out = data.frame(File = character(), Sample = character(),
                          stringsAsFactors = FALSE)

  for (i in seq_len(n.total)) {

    acc   = sra.data$Run[i]
    samp  = sra.data$sample.name[i]
    layout = if ("LibraryLayout" %in% names(sra.data)) {
               sra.data$LibraryLayout[i]
             } else "PAIRED"
    is.paired = toupper(layout) == "PAIRED"

    if (!quiet) message(sprintf("[%d/%d] %s  (%s)", i, n.total, samp, acc))

    # Completion sentinel — fast skip on re-runs.
    # Stored in .sra-sentinels/ so it is never mistaken for a read file.
    sentinel = file.path(sentinels.dir, paste0(samp, "_COMPLETED"))
    if (file.exists(sentinel)) {
      if (!quiet) message("  already completed — skipping")
      rename.out = rbind(rename.out,
                         data.frame(File   = paste0(samp, "_L001"),
                                    Sample = samp, stringsAsFactors = FALSE))
      next
    }

    # Destination paths
    r1.dest = file.path(output.directory, paste0(samp, "_L001_READ1.fastq.gz"))
    r2.dest = if (is.paired)
                file.path(output.directory, paste0(samp, "_L001_READ2.fastq.gz"))
              else NULL

    # Skip if files already exist (prior run without sentinel)
    if (file.exists(r1.dest) && (!is.paired || file.exists(r2.dest))) {
      if (!quiet) message("  output files exist — skipping")
      rename.out = rbind(rename.out,
                         data.frame(File   = paste0(samp, "_L001"),
                                    Sample = samp, stringsAsFactors = FALSE))
      next
    }

    # Build ENA URLs and download
    urls  = .ena.urls(acc, layout)

    r1.ok = .dl(urls["r1"], r1.dest, max.retries, retry.delay, quiet)
    if (!r1.ok) {
      if (file.exists(r1.dest)) file.remove(r1.dest)
      msg = sprintf("  READ1 download failed for %s after %d attempts", acc, max.retries)
      if (!skip.not.found) stop(msg) else { warning(msg); next }
    }

    if (is.paired) {
      r2.ok = .dl(urls["r2"], r2.dest, max.retries, retry.delay, quiet)
      if (!r2.ok) {
        if (file.exists(r1.dest)) file.remove(r1.dest)
        if (file.exists(r2.dest)) file.remove(r2.dest)
        msg = sprintf("  READ2 download failed for %s after %d attempts", acc, max.retries)
        if (!skip.not.found) stop(msg) else { warning(msg); next }
      }
    }

    # Write completion sentinel (in .sra-sentinels/ subdirectory)
    writeLines(c(
      paste0("accession: ", acc),
      paste0("sample:    ", samp),
      paste0("layout:    ", layout),
      paste0("completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    ), sentinel)

    if (!quiet) message("  done")

    rename.out = rbind(rename.out,
                       data.frame(File   = paste0(samp, "_L001"),
                                  Sample = samp, stringsAsFactors = FALSE))
  } # end for i

  # ── Write rename CSV ─────────────────────────────────────────────────────────
  write.csv(rename.out,
            file      = "file_rename_sra.csv",
            row.names = FALSE,
            quote     = FALSE)

  if (!quiet)
    message("\nDone. ", nrow(rename.out), " sample(s) recorded in file_rename_sra.csv.")

  invisible(rename.out)

} # end sraDownload
