# Internal depth helpers for workflow 3.

# Checks the depth-filter settings and returns them in a list.
.validateDepthSettings = function(depth.filter.mode = "site",
                                  min.site.depth = 1,
                                  min.mean.depth = 1,
                                  max.n.proportion = NULL) {

  mode = match.arg(tolower(depth.filter.mode), c("none", "site", "mean", "both"))
  for (name in c("min.site.depth", "min.mean.depth")) {
    value = get(name)
    if (length(value) != 1 || !is.finite(value) || value < 0) {
      stop(name, " must be one non-negative finite value.")
    }
  }
  if (!is.null(max.n.proportion) &&
      (length(max.n.proportion) != 1 || !is.finite(max.n.proportion) ||
       max.n.proportion < 0 || max.n.proportion > 1)) {
    stop("max.n.proportion must be NULL or one value from 0 to 1.")
  }

  list(mode = mode,
       min.site.depth = min.site.depth,
       min.mean.depth = min.mean.depth,
       max.n.proportion = max.n.proportion)
}

#' Calculate full-span per-base depth for workflow 3
#' @param mapping.directory Mapped sample directory.
#' @param output.directory Owned depth-table directory.
#' @param sample.names Samples to scan.
#' @param use.base.recalibration Use recalibrated BAMs.
#' @param samtools.path samtools path/directory or NULL.
#' @param overwrite Recompute cached depth.
#' @param quiet Suppress samtools output while retaining logs.
#' @return Named paths to depth tables.
#' @export
calculateSampleDepth = function(mapping.directory,
                                output.directory,
                                sample.names = NULL,
                                use.base.recalibration = FALSE,
                                samtools.path = NULL,
                                overwrite = FALSE,
                                quiet = TRUE) {

  #Quick checks
  if (!dir.exists(mapping.directory)) { stop("Mapping directory not found.") }
  if (is.null(sample.names)) {
    sample.names = list.dirs(mapping.directory, recursive = FALSE, full.names = FALSE)
  }
  if (length(sample.names) == 0) { stop("No samples are available for depth calculation.") }

  samtools = .toolCommand("samtools", samtools.path)
  .ensureDirectory(output.directory)
  .ensureDirectory("logs/sample_logs")
  paths = setNames(file.path(output.directory, paste0(sample.names, ".depth.tsv")), sample.names)

  for (sample in sample.names) {
    bam = .selectedSampleBam(mapping.directory, sample, use.base.recalibration)
    reference = file.path(mapping.directory, sample, "index", "reference.fa")
    fai = paste0(reference, ".fai")
    if (!file.exists(fai)) { stop("Reference index not found for ", sample) }

    #The signature records the BAM, reference, and samtools flags. A cached depth
    #table is reused only when the signature still matches.
    record = paste0(paths[sample], ".complete")
    info = file.info(bam)
    signature = c(paste0("bam=", normalizePath(bam)),
                  paste0("size=", info$size),
                  paste0("mtime=", as.numeric(info$mtime)),
                  paste0("reference=", normalizePath(reference)),
                  "flags=-aa -s -q 0 -Q 0 -G 0x800")
    reusable = overwrite == FALSE && file.exists(paths[sample]) &&
      file.info(paths[sample])$size > 0 && file.exists(record) &&
      identical(readLines(record, warn = FALSE), signature)
    if (reusable) { next }

    file.remove(c(paths[sample], record))
    log = file.path("logs/sample_logs", sample, "depth.stderr.log")
    .ensureDirectory(dirname(log))
    cmd = paste(samtools, "depth -aa -s -q 0 -Q 0 -G 0x800", shQuote(bam), ">", shQuote(paths[sample]))
    .runCommand(cmd, quiet, "samtools depth", keep.stdout = TRUE, stderr.log = log)
    if (!file.exists(paths[sample]) || file.info(paths[sample])$size == 0) {
      stop("Depth output is empty for ", sample)
    }

    #Validates the shape, contigs, full lengths, and numeric coordinates without
    #loading the BAM.
    depth.table = data.table::fread(paths[sample], header = FALSE, select = 1:3,
                                    col.names = c("contig", "position", "depth"),
                                    showProgress = FALSE)
    index.table = data.table::fread(fai, header = FALSE, select = 1:2,
                                    col.names = c("contig", "length"),
                                    showProgress = FALSE)
    if (any(!depth.table$contig %in% index.table$contig) ||
        any(!is.finite(depth.table$position)) ||
        any(!is.finite(depth.table$depth)) || any(depth.table$depth < 0)) {
      stop("Invalid or truncated depth table for ", sample)
    }

    n.by = table(depth.table$contig)
    counts = data.frame(contig = names(n.by),
                        n = as.integer(n.by),
                        min = as.numeric(tapply(depth.table$position, depth.table$contig, min)),
                        max = as.numeric(tapply(depth.table$position, depth.table$contig, max)))
    check = merge(index.table, counts, by = "contig", all.x = TRUE)
    if (any(is.na(check$n)) || any(check$n != check$length) ||
        any(check$min != 1) || any(check$max != check$length)) {
      stop("Depth table does not cover every indexed reference position for ", sample)
    }
    writeLines(signature, record)
  }#end sample loop

  paths
}

# Restores reference names from full-contig GATK FASTA headers.
.referenceContigNames = function(seqs) {
  headers = sub("^[0-9]+ ", "", names(seqs))
  intervals = regexec(":([0-9]+)-([0-9]+)$", headers)
  parts = regmatches(headers, intervals)
  has.interval = lengths(parts) == 3
  if (any(has.interval)) {
    starts = as.numeric(vapply(parts[has.interval], `[`, character(1), 2))
    ends = as.numeric(vapply(parts[has.interval], `[`, character(1), 3))
    if (any(starts != 1 | ends != Biostrings::width(seqs)[has.interval])) {
      stop("GATK FASTA has a partial contig interval.")
    }
    headers[has.interval] = sub(":1-[0-9]+$", "", headers[has.interval])
  }
  if (anyDuplicated(headers)) {
    stop("GATK FASTA has duplicate contig names.")
  }
  headers
}

# Masks low-depth sites, removes contigs that fail the depth or N rules, and
# writes a per-contig report. Returns the retained sequences.
.filterDepthSequences = function(seqs, depth.file, settings, report.file) {

  depth.table = data.table::fread(depth.file, header = FALSE,
                                  col.names = c("contig", "position", "depth"),
                                  showProgress = FALSE)

  #Use the reference names to match each sequence to its depth records.
  normalized = .referenceContigNames(seqs)
  if (anyDuplicated(normalized) || !setequal(normalized, unique(depth.table$contig))) {
    stop("FASTA and depth contig names do not agree.")
  }
  names(seqs) = normalized

  rows = vector("list", length(seqs))
  keep = rep(TRUE, length(seqs))

  for (i in seq_along(seqs)) {
    contig.depth = depth.table[depth.table$contig == names(seqs)[i]]
    len = Biostrings::width(seqs[i])
    if (nrow(contig.depth) != len) {
      stop("Depth length mismatch for contig ", names(seqs)[i])
    }

    old.n = as.integer(Biostrings::letterFrequency(seqs[i], "N"))
    mean.depth = mean(contig.depth$depth)
    below = sum(contig.depth$depth < settings$min.site.depth)
    newly = 0L

    #Site masking: sites below the threshold become N. An existing N is preserved.
    if (settings$mode %in% c("site", "both")) {
      pos = contig.depth$position[contig.depth$depth < settings$min.site.depth]
      before = as.character(seqs[i])
      if (length(pos)) {
        seqs[i] = Biostrings::replaceAt(seqs[i], IRanges::IRanges(pos, width = 1),
                                        Biostrings::DNAStringSet(rep("N", length(pos))))
      }
      newly = sum(substring(before, pos, pos) != "N")
    }

    final.n = as.integer(Biostrings::letterFrequency(seqs[i], "N"))
    prop = if (len) final.n / len else 1

    #Collects every reason the contig fails, so each removal is recorded once
    reasons = character()
    if (settings$mode %in% c("mean", "both") && mean.depth < settings$min.mean.depth) {
      reasons = c(reasons, "mean_depth")
    }
    filtering = settings$mode != "none" || !is.null(settings$max.n.proportion)
    if (filtering && (len == 0 || final.n == len)) {
      reasons = c(reasons, "no_sequence")
    }
    if (!is.null(settings$max.n.proportion) && prop > settings$max.n.proportion) {
      reasons = c(reasons, "excess_N")
    }

    keep[i] = !length(reasons)
    rows[[i]] = data.frame(contig = names(seqs)[i],
                           length = len,
                           mean_depth = mean.depth,
                           bases_below_threshold = below,
                           preexisting_N = old.n,
                           newly_masked = newly,
                           final_N_proportion = prop,
                           retained = keep[i],
                           reason = paste(unique(reasons), collapse = ";"))
  }#end contig loop

  data.table::fwrite(data.table::rbindlist(rows), report.file, sep = "\t")
  seqs[keep]
}
