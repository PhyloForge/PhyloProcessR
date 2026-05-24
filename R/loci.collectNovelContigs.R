#' @title collectNovelContigs
#'
#' @description Collects and consolidates the per-sample contig FASTAs produced by
#' \code{assembleSharedRegions} into a single to-align FASTA suitable for
#' \code{alignTargets}. Because \code{assembleSharedRegions} names each contig
#' \code{region_contig_N} (e.g. \code{chr3_450000_450800_contig_1}), the
#' locus-to-contig assignment is already encoded in the name — no BLAST against
#' the probe set is needed or used. For each (region, sample) pair the longest
#' assembled contig is selected as the representative sequence. Output contig names
#' follow the \code{locus_|_sample} convention expected by \code{alignTargets}.
#'
#' A per-sample summary CSV is written to \code{logs/novel_contig_summary.csv}.
#'
#' @param contig.directory path to the directory of per-sample contig FASTA files
#' produced by \code{assembleSharedRegions} (or \code{filterHeterozygosity}).
#'
#' @param output.name path prefix for the to-align FASTA. A suffix
#' \code{_to-align.fa} is appended automatically, e.g.
#' \code{data-analysis/novel-loci} writes
#' \code{data-analysis/novel-loci_to-align.fa}.
#'
#' @param min.contig.length minimum contig length in bp. Contigs shorter than
#' this are discarded. Default 200.
#'
#' @param min.taxa minimum number of samples that must have a contig for a locus
#' to be retained. Default 4.
#'
#' @param threads number of CPU threads for parallel sample processing. Default 1.
#'
#' @param overwrite logical. If TRUE, overwrites an existing to-align FASTA.
#' Default FALSE.
#'
#' @return Invisibly returns a data.frame summary. Writes:
#' \itemize{
#'   \item \code{<output.name>_to-align.fa} — combined FASTA for \code{alignTargets}
#'   \item \code{logs/novel_contig_summary.csv} — per-sample contig counts
#' }
#'
#' @export

collectNovelContigs = function(contig.directory = NULL,
                               output.name      = NULL,
                               min.contig.length = 200,
                               min.taxa          = 4,
                               threads           = 1,
                               overwrite         = FALSE) {

  # contig.directory  = "data-analysis/contigs/9_genome-contigs"
  # output.name       = "data-analysis/novel-loci"
  # min.contig.length = 200
  # min.taxa          = 4
  # threads           = 1
  # overwrite         = FALSE

  #Input checks
  if (is.null(contig.directory)) { print("contig.directory not provided."); return(NULL) }
  if (is.null(output.name))      { print("output.name not provided."); return(NULL) }

  out.fa = paste0(output.name, "_to-align.fa")
  if (file.exists(out.fa) && overwrite == FALSE) {
    print(paste0("Output file already exists: ", out.fa,
                 ". Set overwrite = TRUE to redo."))
    return(invisible(NULL))
  }

  dir.create("logs", recursive = TRUE, showWarnings = FALSE)

  fa.files    = list.files(contig.directory, pattern = "\\.fa$", full.names = TRUE)
  sample.names = gsub("\\.fa$", "", basename(fa.files))

  if (length(fa.files) == 0) {
    print("No contig FASTA files found in contig.directory.")
    return(invisible(NULL))
  }

  print(paste0("Collecting novel contigs from ", length(fa.files), " samples..."))

  ##########################################################################
  # Per-sample: read contigs, group by region, keep longest per region
  ##########################################################################
  results = parallel::mclapply(seq_along(fa.files), function(i) {
    tryCatch({

      samp  = sample.names[i]
      ctgs  = Biostrings::readDNAStringSet(fa.files[i])

      if (length(ctgs) == 0) {
        return(data.frame(Sample = samp, InputContigs = 0L,
                          RegionsRepresented = 0L, stringsAsFactors = FALSE))
      }

      # Drop short contigs
      ctgs = ctgs[Biostrings::width(ctgs) >= min.contig.length]

      if (length(ctgs) == 0) {
        return(data.frame(Sample = samp, InputContigs = length(ctgs),
                          RegionsRepresented = 0L, stringsAsFactors = FALSE))
      }

      # Extract region name: everything up to the last "_contig_" occurrence
      # e.g. "chr3_450000_450800_contig_2" → "chr3_450000_450800"
      region.keys = sub("_contig_[0-9]+$", "", names(ctgs))
      # Sanitize to match novel_targets.fa sequence names (handles contigs assembled
      # before scaffold-name sanitization was applied)
      region.keys = gsub("[^A-Za-z0-9_.]", "_", region.keys)

      # For each region keep the longest contig
      best = tapply(seq_along(ctgs), region.keys, function(idx) {
        idx[which.max(Biostrings::width(ctgs)[idx])]
      })
      best.idx = unlist(best, use.names = FALSE)

      out.seqs = ctgs[best.idx]
      # Name as locus_|_sample for alignTargets
      names(out.seqs) = paste0(region.keys[best.idx], "_|_", samp)

      return(list(
        seqs    = out.seqs,
        summary = data.frame(
          Sample             = samp,
          InputContigs       = length(ctgs),
          RegionsRepresented = length(best.idx),
          stringsAsFactors   = FALSE
        )
      ))

    }, error = function(e) {
      warning(sample.names[i], " failed: ", conditionMessage(e))
      return(NULL)
    })
  }, mc.cores = threads)

  ##########################################################################
  # Assemble combined FASTA and summary table
  ##########################################################################
  all.seqs = Biostrings::DNAStringSet()
  sum.rows = list()

  for (r in results) {
    if (is.null(r)) next
    if (is.data.frame(r)) {
      # Sample that had no contigs passing filters
      sum.rows[[length(sum.rows) + 1]] = r
      next
    }
    all.seqs = append(all.seqs, r$seqs)
    sum.rows[[length(sum.rows) + 1]] = r$summary
  }

  summary.df = do.call(rbind, sum.rows)

  if (length(all.seqs) == 0) {
    print("No contigs passed filters. Check min.contig.length and input data.")
    return(invisible(summary.df))
  }

  # Apply min.taxa filter: only keep loci present in >= min.taxa samples
  locus.names = sub("_\\|_.*$", "", names(all.seqs))
  locus.counts = table(locus.names)
  keep.loci = names(locus.counts)[locus.counts >= min.taxa]

  n.before = length(unique(locus.names))
  all.seqs = all.seqs[locus.names %in% keep.loci]
  n.after  = length(keep.loci)

  print(paste0("Loci before min.taxa filter (>= ", min.taxa, " samples): ", n.before))
  print(paste0("Loci retained after filter: ", n.after))
  print(paste0("Total sequences in to-align FASTA: ", length(all.seqs)))

  # Write combined FASTA
  final.loci = as.list(as.character(all.seqs))
  PhyloProcessR::writeFasta(
    sequences = final.loci,
    names     = names(final.loci),
    file.out  = out.fa,
    nbchar    = 1000000,
    as.string = TRUE,
    open      = "w"
  )

  # Per-locus summary: one row per locus with sample count and filter status
  locus.df = data.frame(
    Locus          = names(locus.counts),
    N_samples      = as.integer(locus.counts),
    Passed_filter  = names(locus.counts) %in% keep.loci,
    stringsAsFactors = FALSE
  )
  locus.df = locus.df[order(-locus.df$N_samples), ]

  # Per-sample summary with contig length stats
  # Recompute widths for sequences that passed the min.taxa filter
  if (length(all.seqs) > 0) {
    pass.names  = sub("_\\|_.*$", "", names(all.seqs))
    pass.widths = Biostrings::width(all.seqs)
    per.samp.loci = tapply(pass.names, sub("^.*_\\|_", "", names(all.seqs)), function(x) length(x))
    summary.df$LociPassingFilter = as.integer(per.samp.loci[summary.df$Sample])
    summary.df$LociPassingFilter[is.na(summary.df$LociPassingFilter)] = 0L
  } else {
    summary.df$LociPassingFilter = 0L
  }

  # Write logs
  write.csv(summary.df, file = "logs/novel_contig_summary.csv",  row.names = FALSE)
  write.csv(locus.df,   file = "logs/novel_locus_summary.csv",   row.names = FALSE)

  print(paste0("To-align FASTA written to: ", out.fa))
  print(paste0("Per-sample summary:  logs/novel_contig_summary.csv"))
  print(paste0("Per-locus summary:   logs/novel_locus_summary.csv"))

  return(invisible(summary.df))

}# end function

#END SCRIPT
