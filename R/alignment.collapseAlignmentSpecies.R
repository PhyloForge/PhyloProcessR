#' @title collapseAlignmentSpecies
#'
#' @description Collapses per-individual sequence alignments to species-level
#'   alignments by merging all sequences that share the same \code{Genus_species}
#'   prefix into a single \code{Genus_species_COMBINED} consensus sequence.
#'   The consensus is built column by column: informative (non-gap, non-N)
#'   bases take priority over gaps and ambiguities; when multiple individuals
#'   have conflicting informative bases the first is used.  Collapsed alignments
#'   are written to \code{output.dir} in phylip format.
#'
#'   Accepts either a directory of pre-existing phylip / fasta alignments
#'   (\code{alignment.dir}) or a concatenated NEXUS file (\code{nexus.file}).
#'   When a NEXUS file is supplied it is first split into per-locus phylip files
#'   by \code{convertNexusPartitions}; those files then become the input.
#'
#' @param alignment.dir path to a directory of alignment files to collapse.
#'   Ignored when \code{nexus.file} is supplied and \code{nexus.output.dir} is
#'   set, in which case the per-locus phylip files written there are used.
#'
#' @param alignment.format format of files in \code{alignment.dir}: \code{"phylip"}
#'   (default) or \code{"fasta"}.
#'
#' @param output.dir path to the directory where collapsed alignments are saved.
#'   Default \code{"constrained-alignments"}.
#'
#' @param nexus.file optional path to a concatenated NEXUS file with a
#'   \code{BEGIN SETS / charset} block.  When supplied, the NEXUS is split into
#'   per-locus phylip files before collapsing.
#'
#' @param nexus.output.dir directory where the per-locus phylip files produced
#'   from \code{nexus.file} are written (and subsequently read).  Required when
#'   \code{nexus.file} is supplied.
#'
#' @param max.missing.percent maximum percentage of gap / N / ? characters
#'   allowed in a collapsed (merged) species sequence.  Species whose merged
#'   sequence exceeds this threshold are dropped from that alignment.  Default
#'   100 (keep all).
#'
#' @param min.taxa.alignment minimum number of collapsed species required to
#'   save an alignment.  Default 4.
#'
#' @param overwrite if \code{TRUE} delete and recreate \code{output.dir} before
#'   writing.  Default \code{FALSE}.
#'
#' @param threads number of parallel processing threads.  Default 1.
#'
#' @param memory total memory in GB to allocate across threads.  Default 1.
#'
#' @param quiet if \code{TRUE} suppress per-alignment progress messages.
#'   Default \code{TRUE}.
#'
#' @return invisibly saves collapsed phylip alignments to \code{output.dir};
#'   nothing is returned to R.
#'
#' @export

collapseAlignmentSpecies = function(alignment.dir       = NULL,
                                    alignment.format     = "phylip",
                                    output.dir           = "constrained-alignments",
                                    nexus.file           = NULL,
                                    nexus.output.dir     = NULL,
                                    max.missing.percent  = 100,
                                    min.taxa.alignment   = 4,
                                    overwrite            = FALSE,
                                    threads              = 1,
                                    memory               = 1,
                                    quiet                = TRUE) {

  # Debug
  # alignment.dir      = "data-analysis/legacy-integration/trimmed_legacy-only"
  # alignment.format   = "phylip"
  # output.dir         = "constrained-alignments"
  # nexus.file         = "legacy/Glassfrog_matrix_2026.nex"
  # nexus.output.dir   = "legacy/legacy-alignments"
  # max.missing.percent = 100
  # min.taxa.alignment  = 4

  ##############################################################################
  # 0. Parameter checks
  ##############################################################################
  if (is.null(nexus.file) && is.null(alignment.dir)) {
    stop("Error: supply either alignment.dir or nexus.file (or both).")
  }
  if (!is.null(nexus.file) && is.null(nexus.output.dir)) {
    stop("Error: nexus.output.dir is required when nexus.file is supplied.")
  }

  ##############################################################################
  # 1. NEXUS → per-locus phylip conversion (optional)
  ##############################################################################
  if (!is.null(nexus.file)) {
    if (!file.exists(nexus.file)) {
      stop(paste0("Error: nexus.file not found: ", nexus.file))
    }
    print(paste0("Splitting NEXUS → per-locus phylip in: ", nexus.output.dir))
    PhyloProcessR::convertNexusPartitions(
      nexus.file          = nexus.file,
      output.directory    = nexus.output.dir,
      output.format       = "phylip",
      max.missing.percent = max.missing.percent,
      overwrite           = overwrite,
      quiet               = quiet
    )
    # Use the converted files as the alignment source
    alignment.dir    = nexus.output.dir
    alignment.format = "phylip"
  }

  ##############################################################################
  # 2. Set up output directory
  ##############################################################################
  if (alignment.dir == output.dir) {
    stop("Error: output.dir must differ from alignment.dir.")
  }
  if (dir.exists(output.dir)) {
    if (overwrite) {
      unlink(output.dir, recursive = TRUE)
      dir.create(output.dir, recursive = TRUE)
    }
    # else: silently reuse; individual files are skipped if already present
  } else {
    dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  }

  ##############################################################################
  # 3. List alignment files
  ##############################################################################
  ext.pattern = if (alignment.format == "phylip") "\\.phy(lip)?$" else "\\.fa(sta)?$"
  align.files = list.files(alignment.dir, pattern = ext.pattern,
                           full.names = FALSE, recursive = FALSE)
  if (length(align.files) == 0) {
    stop(paste0("Error: no alignment files found in: ", alignment.dir))
  }
  print(paste0("collapseAlignmentSpecies: processing ", length(align.files),
               " alignment(s) → ", output.dir))

  ##############################################################################
  # 4. Helper: species key from a tip label
  #    First two underscore-delimited tokens = Genus_species.
  #    If the second token is cf/aff/sp, include three tokens.
  ##############################################################################
  sp.key = function(x) {
    parts = strsplit(x, "_")[[1L]]
    n = if (length(parts) >= 3L && parts[2L] %in% c("cf", "aff", "sp"))
      min(3L, length(parts))
    else
      min(2L, length(parts))
    paste(parts[seq_len(n)], collapse = "_")
  }

  ##############################################################################
  # 5. Helper: column-by-column consensus of N sequences
  ##############################################################################
  col.consensus = function(seqs.char) {
    # seqs.char: character vector of equal-length strings
    mat = matrix(unlist(strsplit(seqs.char, "")),
                 nrow = length(seqs.char), byrow = TRUE)
    apply(mat, 2L, function(col) {
      if (length(unique(col)) == 1L) return(col[1L])
      informative = col[!col %in% c("-", "N", "n", "?")]
      if (length(informative) == 0L) "-" else informative[1L]
    })
  }

  ##############################################################################
  # 6. Per-alignment collapse loop (parallel)
  ##############################################################################
  parallel::mclapply(seq_along(align.files), function(i) {
  tryCatch({

    aln.path  = file.path(alignment.dir, align.files[i])
    save.name = sub("\\.(phy(lip)?|fa(sta)?)$", "", align.files[i])
    out.path  = file.path(output.dir, paste0(save.name, ".phy"))

    # Skip already-done files when overwrite = FALSE
    if (!overwrite && file.exists(out.path)) {
      if (!quiet) print(paste0(save.name, " already exists, skipping."))
      return(NULL)
    }

    # Load alignment
    aln = if (alignment.format == "phylip") {
      tryCatch({
        raw = Biostrings::readDNAMultipleAlignment(aln.path, format = "phylip")
        Biostrings::DNAStringSet(raw)
      }, error = function(e) {
        # Some phylip files written with '?' — sanitise and retry
        tmp = tempfile(fileext = ".phy")
        writeLines(gsub("\\?", "N", readLines(aln.path, warn = FALSE)), tmp)
        on.exit(unlink(tmp), add = TRUE)
        raw = Biostrings::readDNAMultipleAlignment(tmp, format = "phylip")
        Biostrings::DNAStringSet(raw)
      })
    } else {
      Biostrings::readDNAStringSet(aln.path)
    }

    if (length(aln) == 0) return(NULL)

    tip.nms  = names(aln)
    sp.keys  = sapply(tip.nms, sp.key, USE.NAMES = FALSE)
    unique.sp = unique(sp.keys)

    # Collapse each species
    collapsed = Biostrings::DNAStringSet()

    for (sp in unique.sp) {
      idx  = which(sp.keys == sp)
      seqs = as.character(aln[idx])

      if (length(seqs) == 1L) {
        # Single individual — just rename
        consensus.str = seqs
      } else {
        # Multiple individuals — column-by-column merge
        consensus.str = paste(col.consensus(seqs), collapse = "")
      }

      # Apply max.missing.percent filter on the merged sequence
      n.total = nchar(consensus.str)
      n.gap   = nchar(gsub("[^-?nN]", "", consensus.str, ignore.case = TRUE))
      pct.gap = if (n.total > 0) 100 * n.gap / n.total else 100
      if (pct.gap > max.missing.percent) next

      sp.seq = Biostrings::DNAStringSet(consensus.str)
      names(sp.seq) = paste0(sp, "_COMBINED")
      collapsed = append(collapsed, sp.seq)
    }

    # Minimum taxa filter
    if (length(collapsed) < min.taxa.alignment) {
      if (!quiet) {
        print(paste0(save.name, ": only ", length(collapsed),
                     " species after collapse (< min.taxa.alignment = ",
                     min.taxa.alignment, ") — skipping."))
      }
      return(NULL)
    }

    # Write phylip
    write.temp  = strsplit(as.character(collapsed), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp))
    PhyloProcessR::writePhylip(alignment  = aligned.set,
                               file       = out.path,
                               interleave = FALSE,
                               strict     = FALSE)

    if (!quiet) {
      print(paste0(save.name, ": ", length(aln), " individuals → ",
                   length(collapsed), " species (_COMBINED)"))
    }

    return(invisible(NULL))

  }, error = function(e) {
    warning(align.files[i], " failed: ", conditionMessage(e))
    NULL
  })
  }, mc.cores = threads) #end mclapply

  print(paste0("collapseAlignmentSpecies: done. Output in: ", output.dir))
  return(invisible(NULL))

}#end function
