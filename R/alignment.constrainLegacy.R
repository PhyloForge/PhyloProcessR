#' @title constrainLegacy
#'
#' @description Creates an IQ-TREE backbone constraint tree from a sequence
#'   capture phylogeny.  Two modes are available:
#'
#'   \strong{Default mode} (\code{collapse.to.species = FALSE}): fuzzy-matches
#'   capture tip labels to the taxon names used in an integrated legacy+capture
#'   dataset and renames tips accordingly.
#'
#'   \strong{Species-collapse mode} (\code{collapse.to.species = TRUE}): strips
#'   the voucher ID from every capture tip to obtain the species name, selects
#'   the best-sequenced individual per species, prunes all other conspecific
#'   tips, and renames the surviving tip to \code{Genus_species_COMBINED}.  The
#'   resulting tree can be used as a backbone constraint for a species-level
#'   supermatrix analysis.  If \code{alignment.dir} is supplied, the best
#'   individual is the one with the most total non-gap bases across all
#'   alignment files; otherwise the individual with the shortest terminal branch
#'   (a proxy for sequence completeness) is chosen.
#'
#' @param capture.tree a \code{phylo} object or a path to a Newick / NEXUS file
#'   containing the sequence-capture phylogeny whose topology will be used as
#'   the constraint.
#'
#' @param legacy.taxa character vector of taxon names as they appear in the
#'   integrated legacy+capture alignment.  A \code{phylo} object is also
#'   accepted (its \code{tip.label} slot is used).  In species-collapse mode
#'   this is used only to verify that each collapsed species is represented in
#'   the dataset; it may be \code{NULL} to skip the check.
#'
#' @param output.file path for the output constraint tree (Newick format).  If
#'   \code{NULL} the tree is returned but not written.
#'
#' @param name.match character; matching strategy for the default (non-collapse)
#'   mode.
#'   \describe{
#'     \item{"fuzzy"}{Strip all separators (\code{[-_. ]}) and lowercase
#'       (default).}
#'     \item{"species"}{Remove the last underscore-delimited token (voucher ID)
#'       and match on \code{Genus_species} only.}
#'     \item{"exact"}{Require identical strings.}
#'   }
#'
#' @param collapse.to.species logical; if \code{TRUE} activate species-collapse
#'   mode (see Description).  Default \code{FALSE}.
#'
#' @param alignment.dir path to a directory of alignment files.  Used in
#'   species-collapse mode to select the individual with the most non-gap bases
#'   across all loci.  If \code{NULL} the individual with the shortest terminal
#'   branch length is chosen instead.
#'
#' @param alignment.format format of alignment files in \code{alignment.dir};
#'   \code{"phylip"} (default) or \code{"fasta"}.
#'
#' @param drop.unmatched logical; if \code{FALSE} (default) the function stops
#'   with an error when any tip cannot be matched.  If \code{TRUE} unmatched
#'   tips are pruned with a warning.
#'
#' @param quiet logical; if \code{TRUE} progress messages are suppressed.
#'
#' @return a \code{phylo} object written to \code{output.file} in Newick format
#'   (returned invisibly).
#'
#' @export

constrainLegacy = function(capture.tree     = NULL,
                           legacy.taxa      = NULL,
                           output.file      = NULL,
                           name.match       = c("fuzzy", "species", "exact"),
                           collapse.to.species = FALSE,
                           alignment.dir    = NULL,
                           alignment.format = "phylip",
                           drop.unmatched   = FALSE,
                           quiet            = FALSE) {

  # Debug
  # capture.tree        = "concat_tree.treefile"
  # legacy.taxa         = names(ape::read.nexus.data("legacy.nex"))
  # output.file         = "constraint.tre"
  # collapse.to.species = TRUE
  # alignment.dir       = "data-analysis/alignments/trimmed_all-markers"

  ##############################################################################
  # 0. Parameter checks
  ##############################################################################
  if (is.null(capture.tree)) { stop("Error: capture.tree is required.") }
  if (!collapse.to.species && is.null(legacy.taxa)) {
    stop("Error: legacy.taxa is required when collapse.to.species = FALSE.")
  }

  name.match = match.arg(name.match)

  ##############################################################################
  # 1. Load capture tree
  ##############################################################################
  if (inherits(capture.tree, "phylo")) {
    tree = capture.tree
  } else if (is.character(capture.tree) && length(capture.tree) == 1) {
    if (!file.exists(capture.tree)) {
      stop(paste0("Error: capture tree file not found: ", capture.tree))
    }
    tree = tryCatch(
      ape::read.tree(capture.tree),
      error = function(e) {
        tryCatch(ape::read.nexus(capture.tree),
                 error = function(e2)
                   stop("Error: could not parse capture.tree as Newick or NEXUS."))
      }
    )
  } else {
    stop("Error: capture.tree must be a phylo object or a file path string.")
  }

  cap.tips = tree$tip.label
  if (length(cap.tips) == 0) { stop("Error: capture tree has no tip labels.") }

  ##############################################################################
  # BRANCH: species-collapse mode
  ##############################################################################
  if (collapse.to.species) {

    # --- 2a. Derive species name for every tip --------------------------------
    # Take the first two underscore-delimited tokens as Genus_species.
    # Exception: if the second token is "cf", "aff", or "sp", include the third
    # token too (e.g. Centrolene_cf_venezuelense_MAR371 → Centrolene_cf_venezuelense).
    # This handles voucher IDs that themselves contain underscores
    # (e.g. MNCN-ADN_51750) which would leave fragments if we simply stripped
    # the last token.
    sp.of.tip = sapply(cap.tips, function(x) {
      parts = strsplit(x, "_")[[1]]
      if (length(parts) >= 3 && parts[2] %in% c("cf", "aff", "sp")) {
        paste(parts[1:min(3, length(parts))], collapse = "_")
      } else {
        paste(parts[1:min(2, length(parts))], collapse = "_")
      }
    }, USE.NAMES = FALSE)
    unique.sp  = unique(sp.of.tip)

    if (!quiet) {
      print(paste0("constrainLegacy (species-collapse): ",
                   length(cap.tips), " tips → ", length(unique.sp), " species"))
    }

    # --- 2b. Tally non-gap bases per specimen (if alignment.dir supplied) ----
    specimen.bases = NULL
    if (!is.null(alignment.dir)) {
      if (!quiet) { print("Counting non-gap bases per specimen across alignments...") }
      ext = if (alignment.format == "phylip") "\\.phy$" else "\\.fa(sta)?$"
      aln.files = list.files(alignment.dir, pattern = ext,
                             full.names = TRUE, recursive = FALSE)
      if (length(aln.files) == 0) {
        warning("No alignment files found in alignment.dir; ",
                "falling back to branch-length selection.")
      } else {
        specimen.bases = setNames(numeric(length(cap.tips)), cap.tips)
        for (af in aln.files) {
          aln = tryCatch({
            if (alignment.format == "phylip") {
              Biostrings::DNAStringSet(
                Biostrings::readDNAMultipleAlignment(af, format = "phylip"))
            } else {
              Biostrings::readDNAStringSet(af)
            }
          }, error = function(e) NULL)
          if (is.null(aln)) next
          for (tip in cap.tips) {
            hit = which(names(aln) == tip)
            if (length(hit) == 0) next
            n.bases = nchar(gsub("[-?nN]", "", as.character(aln[hit[1]]),
                                 ignore.case = TRUE))
            specimen.bases[tip] = specimen.bases[tip] + n.bases
          }
        }
      }
    }

    # Pre-compute terminal branch lengths (fallback selector)
    term.bl = setNames(rep(NA_real_, length(cap.tips)), cap.tips)
    for (ti in seq_along(cap.tips)) {
      edge.row = which(tree$edge[, 2] == ti)
      if (length(edge.row) == 1L)
        term.bl[cap.tips[ti]] = tree$edge.length[edge.row]
    }

    # --- 2c. Pick one representative per species -----------------------------
    tips.to.drop = character(0)
    chosen.tips  = character(0)

    for (sp in unique.sp) {
      sp.tips = cap.tips[sp.of.tip == sp]

      if (length(sp.tips) == 1L) {
        chosen.tips = c(chosen.tips, sp.tips)
        next
      }

      # Multiple individuals: rank by total non-gap bases, or branch length
      if (!is.null(specimen.bases)) {
        best = sp.tips[which.max(specimen.bases[sp.tips])]
      } else {
        # Shortest terminal branch = closest to consensus = proxy for completeness
        bl.vals = term.bl[sp.tips]
        best    = sp.tips[which.min(bl.vals)]
      }

      chosen.tips  = c(chosen.tips, best)
      drop.these   = sp.tips[sp.tips != best]
      tips.to.drop = c(tips.to.drop, drop.these)

      if (!quiet) {
        criterion = if (!is.null(specimen.bases))
          paste0(specimen.bases[best], " bp")
        else
          paste0("branch length ", round(term.bl[best], 5))
        print(paste0(sp, ": keeping '", best, "' (",
                     length(sp.tips), " individuals; ", criterion, ")"))
      }
    }#end for sp

    # --- 2d. Prune duplicate tips and rename to _COMBINED --------------------
    if (length(tips.to.drop) > 0) {
      tree = ape::drop.tip(tree, tips.to.drop)
    }
    # Derive species key for each surviving tip using the same 2-token rule
    tree$tip.label = sapply(tree$tip.label, function(x) {
      parts = strsplit(x, "_")[[1]]
      sp = if (length(parts) >= 3 && parts[2] %in% c("cf", "aff", "sp"))
        paste(parts[1:min(3, length(parts))], collapse = "_")
      else
        paste(parts[1:min(2, length(parts))], collapse = "_")
      paste0(sp, "_COMBINED")
    }, USE.NAMES = FALSE)

    # --- 2e. Optional: verify each species exists in legacy.taxa -------------
    if (!is.null(legacy.taxa)) {
      if (inherits(legacy.taxa, "phylo")) { leg.names = legacy.taxa$tip.label
      } else { leg.names = legacy.taxa }

      # Species-level keys from legacy names
      leg.sp = unique(gsub("_[^_]+$", "", leg.names))
      tree.sp.raw = gsub("_COMBINED$", "", tree$tip.label)
      unmatched.sp = tree.sp.raw[!tree.sp.raw %in% leg.sp]

      if (length(unmatched.sp) > 0) {
        msg = paste0(length(unmatched.sp),
                     " species in constraint tree not found in legacy.taxa:\n",
                     paste0("  ", unmatched.sp, collapse = "\n"))
        if (drop.unmatched) {
          warning(msg)
          tree = ape::drop.tip(tree, paste0(unmatched.sp, "_COMBINED"))
        } else {
          stop(paste0("Error: ", msg,
                      "\nSet drop.unmatched = TRUE to prune them instead."))
        }
      }
    }

    # --- 2f. Write and return -------------------------------------------------
    if (!is.null(output.file)) {
      ape::write.tree(tree, file = output.file)
      if (!quiet) {
        print(paste0("Constraint tree written to: ", output.file,
                     "  (", length(tree$tip.label), " species)"))
      }
    }
    return(invisible(tree))

  }#end collapse.to.species

  ##############################################################################
  # BRANCH: default mode — match capture tips to legacy names
  ##############################################################################

  # 2. Resolve legacy.taxa
  if (inherits(legacy.taxa, "phylo")) {
    leg.names = legacy.taxa$tip.label
  } else {
    leg.names = legacy.taxa
  }
  if (length(leg.names) == 0) { stop("Error: legacy.taxa contains no names.") }

  # 3. Build key function
  if (name.match == "species") {
    key.fn = function(x) gsub("_[^_]+$", "", x)
  } else if (name.match == "fuzzy") {
    key.fn = function(x) tolower(gsub("[-_. ]", "", x))
  } else {
    key.fn = identity
  }

  cap.keys = key.fn(cap.tips)
  leg.keys = key.fn(leg.names)

  # 4. Match each capture tip to exactly one legacy name
  matched.legacy.names = character(length(cap.tips))
  unmatched.tips       = character(0)
  ambiguous.tips       = character(0)

  for (i in seq_along(cap.tips)) {
    hit = which(leg.keys == cap.keys[i])

    if (length(hit) == 1L) {
      matched.legacy.names[i] = leg.names[hit]

    } else if (length(hit) == 0L) {
      matched.legacy.names[i] = NA_character_
      unmatched.tips = c(unmatched.tips, cap.tips[i])

    } else {
      # Multiple hits: prefer exact string match; otherwise take first
      exact.hit = which(leg.names[hit] == cap.tips[i])
      if (length(exact.hit) == 1L) {
        matched.legacy.names[i] = leg.names[hit[exact.hit]]
      } else {
        matched.legacy.names[i] = leg.names[hit[1]]
        ambiguous.tips = c(ambiguous.tips,
                           paste0("  '", cap.tips[i], "' → ",
                                  paste0("'", leg.names[hit], "'", collapse = ", ")))
      }
    }
  }

  # 5. Warn on ambiguous matches
  if (length(ambiguous.tips) > 0 && !quiet) {
    warning(paste0(length(ambiguous.tips),
                   " capture tip(s) matched multiple legacy names; first used:\n",
                   paste(ambiguous.tips, collapse = "\n")))
  }

  # 6. Handle unmatched tips
  if (length(unmatched.tips) > 0) {
    msg = paste0(length(unmatched.tips),
                 " capture tip(s) could not be matched (name.match = '",
                 name.match, "'):\n",
                 paste0("  ", unmatched.tips, collapse = "\n"), "\n",
                 "Set drop.unmatched = TRUE to prune instead of stopping.")
    if (drop.unmatched) {
      warning(paste0("Dropping ", length(unmatched.tips), " unmatched tip(s):\n",
                     paste0("  ", unmatched.tips, collapse = "\n")))
      tree = ape::drop.tip(tree, unmatched.tips)
      keep.idx             = !cap.tips %in% unmatched.tips
      matched.legacy.names = matched.legacy.names[keep.idx]
      cap.tips             = cap.tips[keep.idx]
    } else {
      stop(paste0("Error: ", msg))
    }
  }

  # 7. Guard against duplicate legacy names (invalid tree)
  dup.mapped = duplicated(matched.legacy.names)
  if (any(dup.mapped)) {
    dup.nms    = unique(matched.legacy.names[dup.mapped])
    dup.detail = sapply(dup.nms, function(nm) {
      src = cap.tips[matched.legacy.names == nm]
      paste0("  '", nm, "' ← ", paste0("'", src, "'", collapse = ", "))
    })
    stop(paste0("Error: ", length(dup.nms), " legacy name(s) matched by >1 capture tip:\n",
                paste(dup.detail, collapse = "\n"), "\n",
                "Use a stricter name.match strategy or rename the conflicting taxa."))
  }

  # 8. Rename tips and write
  tree$tip.label = matched.legacy.names

  if (!is.null(output.file)) {
    ape::write.tree(tree, file = output.file)
    if (!quiet) { print(paste0("Constraint tree written to: ", output.file)) }
  }

  if (!quiet) {
    print(paste0("constrainLegacy: ", length(cap.tips), " capture tips → ",
                 length(matched.legacy.names), " legacy names matched",
                 if (length(unmatched.tips) > 0)
                   paste0(" (", length(unmatched.tips), " dropped)")
                 else ""))
  }

  return(invisible(tree))

}#end function
