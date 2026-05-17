#' @title filterHeterozygosity
#'
#' @description Filters per-sample FASTA contig files based on the proportion
#'   of IUPAC ambiguity codes (heterozygous sites). For each contig the
#'   proportion of ambiguous bases (R, Y, K, M, S, W, B, D, H, V) is computed;
#'   contigs at or above the threshold are written to removed.directory while
#'   contigs below the threshold are written to output.directory. Filtering is
#'   per-contig, not per-sample, so a sample with one problematic locus keeps
#'   all its other loci.
#'
#'   High IUPAC density is a good proxy for chimeric assembly or mis-assembled
#'   paralogous copies: when two divergent sequences are merged into one contig
#'   the sites where they differ all become ambiguous. Genuine diploid
#'   heterozygosity is also flagged, so choose the threshold with your taxon in
#'   mind — 5 \% is reasonable for vertebrates but may be too strict for
#'   high-diversity invertebrate groups.
#'
#'   Note: this filter does NOT detect clean contamination from another organism
#'   (a foreign-species contig has zero IUPAC codes and passes). For that,
#'   use removeContamination() on reads before assembly.
#'
#'   Two levels of logging are written:
#'   \itemize{
#'     \item \code{logs/sample_logs/<Sample>_heterozygosity.csv} — one row per
#'       contig with its length, IUPAC count, proportion, exempt flag, and
#'       pass/fail result.
#'     \item \code{logs/filterHeterozygosity_summary.csv} — one row per sample
#'       summarising total, kept, removed, exempt contig counts, percentage
#'       removed, and mean / max IUPAC proportion.
#'   }
#'
#' @param iupac.directory path to the directory of per-sample FASTA files
#'   containing IUPAC-coded sequences (e.g. output of VCFtoContigs() with
#'   \code{ambiguity.codes = TRUE}).
#'
#' @param output.directory path where per-sample FASTA files with only the
#'   passing contigs will be saved.
#'
#' @param removed.directory path where per-sample FASTA files containing only
#'   the removed contigs will be saved (for inspection).
#'
#' @param threshold numeric (0–1); maximum allowed proportion of IUPAC
#'   ambiguity characters per contig. Contigs at or above this value are
#'   flagged. Default: \code{0.05} (5 \%).
#'
#' @param min.length integer; contigs shorter than this value (bp) are exempt
#'   from the proportion filter — their pass/fail is recorded in the per-sample
#'   log but they are always kept. Default: \code{100}.
#'
#' @param threads number of parallel samples to process simultaneously.
#'
#' @param memory total RAM in GB; divided equally across threads.
#'
#' @param overwrite logical; if TRUE existing output and removed directories are
#'   deleted and recreated before processing.
#'
#' @return invisibly; the per-sample summary data.frame is also returned.
#'
#' @export

filterHeterozygosity = function(iupac.directory = NULL,
                                output.directory = NULL,
                                removed.directory = NULL,
                                threshold = 0.05,
                                min.length = 100,
                                threads = 1,
                                memory = 1,
                                overwrite = FALSE) {

  # setwd("/Volumes/LaCie/Mantellidae")
  # iupac.directory   = "data-analysis/contigs/5_iupac-contigs"
  # removed.directory = "data-analysis/contigs/6_removed-contigs"
  # output.directory  = "data-analysis/contigs/7_filtered-contigs"
  # threshold  = 0.05
  # min.length = 100
  # threads    = 4
  # memory     = 20
  # overwrite  = FALSE

  #################################################
  ### Checks and directory setup
  #################################################

  if (is.null(iupac.directory))  { stop("Please provide the iupac.directory.") }
  if (is.null(output.directory)) { stop("Please provide an output.directory.") }
  if (is.null(removed.directory)){ stop("Please provide a removed.directory.") }
  if (iupac.directory == output.directory)  { stop("iupac.directory and output.directory must differ.") }
  if (iupac.directory == removed.directory) { stop("iupac.directory and removed.directory must differ.") }

  if (dir.exists(output.directory)) {
    if (overwrite) { system(paste0("rm -r ", output.directory)); dir.create(output.directory) }
  } else { dir.create(output.directory) }

  if (dir.exists(removed.directory)) {
    if (overwrite) { system(paste0("rm -r ", removed.directory)); dir.create(removed.directory) }
  } else { dir.create(removed.directory) }

  if (!dir.exists("logs/sample_logs")) {
    dir.create("logs/sample_logs", recursive = TRUE, showWarnings = FALSE)
  }

  file.names = list.files(iupac.directory)
  if (length(file.names) == 0) { stop("No FASTA files found in iupac.directory.") }

  #################################################
  ### Per-sample filtering (parallel)
  #################################################

  results = parallel::mclapply(seq_along(file.names), function(i) {
    tryCatch({

      sample.name = gsub("\\.fa$", "", file.names[i])

      contigs = Biostrings::readDNAStringSet(
        paste0(iupac.directory, "/", file.names[i]), format = "fasta"
      )

      if (length(contigs) == 0) {
        warning(file.names[i], ": empty FASTA — skipping.")
        return(data.frame(Sample = sample.name, TotalContigs = 0L,
                          KeptContigs = 0L, RemovedContigs = 0L,
                          PctRemoved = NA_real_, MeanIUPACprop = NA_real_,
                          MaxIUPACprop = NA_real_, MeanContigLength = NA_real_,
                          stringsAsFactors = FALSE))
      }

      widths = Biostrings::width(contigs)
      seqs   = as.character(contigs)

      # Count IUPAC ambiguity bases with a single gsub pass (faster than
      # running str_count 10 times)
      iupac.counts = nchar(gsub("[^RYKMSWBDHV]", "", seqs))
      iupac.prop   = iupac.counts / widths

      # Exempt flag: contig is below min.length — always kept regardless of proportion
      exempt  = widths < min.length

      # Flagged: long enough AND proportion meets/exceeds threshold
      flagged = !exempt & (iupac.prop >= threshold)

      low.contigs  = contigs[!flagged]
      high.contigs = contigs[flagged]

      #------------------------------------------------------
      # Per-contig log for this sample
      #------------------------------------------------------
      contig.log = data.frame(
        Contig     = names(contigs),
        Length     = widths,
        IUPACcount = iupac.counts,
        IUPACprop  = round(iupac.prop, 6),
        Failed     = flagged,
        stringsAsFactors = FALSE
      )
      write.csv(contig.log,
                file = paste0("logs/sample_logs/", sample.name, "_heterozygosity.csv"),
                row.names = FALSE)

      #------------------------------------------------------
      # Write passing and removed FASTA files (only if non-empty
      # to avoid empty files crashing downstream functions)
      #------------------------------------------------------
      if (length(low.contigs) > 0) {
        final.loci = as.list(as.character(low.contigs))
        PhyloProcessR::writeFasta(
          sequences = final.loci, names = names(final.loci),
          paste0(output.directory, "/", file.names[i]),
          nbchar = 1000000, as.string = TRUE
        )
      }

      if (length(high.contigs) > 0) {
        final.loci = as.list(as.character(high.contigs))
        PhyloProcessR::writeFasta(
          sequences = final.loci, names = names(final.loci),
          paste0(removed.directory, "/", file.names[i]),
          nbchar = 1000000, as.string = TRUE
        )
      }

      #------------------------------------------------------
      # Return per-sample summary row
      #------------------------------------------------------
      data.frame(
        Sample           = sample.name,
        TotalContigs     = length(contigs),
        KeptContigs      = length(low.contigs),
        RemovedContigs   = length(high.contigs),
        PctRemoved       = round(length(high.contigs) / length(contigs) * 100, 2),
        MeanIUPACprop    = round(mean(iupac.prop), 6),
        MaxIUPACprop     = round(max(iupac.prop), 6),
        MeanContigLength = round(mean(widths), 1),
        stringsAsFactors = FALSE
      )

    }, error = function(e) {
      warning(file.names[i], " failed: ", conditionMessage(e))
      data.frame(Sample = gsub("\\.fa$", "", file.names[i]),
                 TotalContigs = NA_integer_, KeptContigs = NA_integer_,
                 RemovedContigs = NA_integer_,
                 PctRemoved = NA_real_, MeanIUPACprop = NA_real_,
                 MaxIUPACprop = NA_real_, MeanContigLength = NA_real_,
                 stringsAsFactors = FALSE)
    })
  }, mc.cores = threads)

  #################################################
  ### Write cross-sample summary CSV
  ### (done after the parallel loop to avoid race conditions)
  #################################################

  summary.df = do.call(rbind, results[!sapply(results, is.null)])

  if (!is.null(summary.df) && nrow(summary.df) > 0) {

    # Append-not-overwrite pattern: merge with any existing summary so that
    # partial re-runs (overwrite = FALSE) accumulate correctly
    out.csv = "logs/filterHeterozygosity_summary.csv"
    if (file.exists(out.csv)) {
      existing = read.csv(out.csv, stringsAsFactors = FALSE)
      existing = existing[!existing$Sample %in% summary.df$Sample, ]
      summary.df = rbind(existing, summary.df)
    }
    write.csv(summary.df, file = out.csv, row.names = FALSE)

    # Print totals
    good.rows = summary.df[!is.na(summary.df$TotalContigs), ]
    cat(" filterHeterozygosity complete:\n",
        " ", sum(good.rows$KeptContigs),    "contigs kept\n",
        " ", sum(good.rows$RemovedContigs), "contigs removed\n",
        "  across", nrow(good.rows), "samples.\n",
        "  Summary: logs/filterHeterozygosity_summary.csv\n",
        "  Per-contig detail: logs/sample_logs/<Sample>_heterozygosity.csv\n")
  }

  invisible(summary.df)

} # end function
