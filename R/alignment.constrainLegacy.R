#' @title constrainLegacy
#'
#' @description Creates an IQ-TREE backbone constraint tree from a sequence
#'   capture phylogeny by fuzzy-matching capture tip labels to the taxon names
#'   used in an integrated legacy+capture dataset. The capture tree topology is
#'   preserved; only tip labels are replaced with their matched legacy names so
#'   the constraint tree is compatible with the final supermatrix alignment.
#'
#'   By default the function stops with an informative error if any capture tip
#'   cannot be matched to a legacy name, forcing the user to resolve naming
#'   discrepancies before proceeding.  Set \code{drop.unmatched = TRUE} to
#'   silently drop unresolvable tips instead.
#'
#' @param capture.tree a \code{phylo} object or a path to a Newick / NEXUS file
#'   containing the sequence-capture phylogeny whose topology will be used as
#'   the constraint.
#'
#' @param legacy.taxa character vector of taxon names as they appear in the
#'   integrated legacy+capture alignment (i.e. the tip labels that IQ-TREE will
#'   see).  A \code{phylo} object is also accepted; its \code{tip.label} slot is
#'   used.
#'
#' @param output.file path for the output constraint tree (written in Newick
#'   format).  If NULL the tree is returned but not written.
#'
#' @param name.match character; matching strategy used to link capture names to
#'   legacy names.
#'   \describe{
#'     \item{"fuzzy"}{Strip all separators (\code{[-_. ]}) and lowercase before
#'       comparing.  Handles minor punctuation and capitalisation differences
#'       between datasets (default).}
#'     \item{"species"}{Remove the last underscore-delimited token (voucher ID)
#'       from each name and match on \code{Genus_species} only.}
#'     \item{"exact"}{Require identical strings.}
#'   }
#'
#' @param drop.unmatched logical; if \code{FALSE} (default) the function stops
#'   with an error listing every capture tip that could not be matched, so the
#'   user can fix naming issues.  If \code{TRUE} unmatched tips are pruned from
#'   the constraint tree and a warning is issued.
#'
#' @param quiet logical; if \code{TRUE} progress messages are suppressed.
#'
#' @return a \code{phylo} object with tip labels replaced by matched legacy
#'   names, written to \code{output.file} in Newick format (invisibly).
#'
#' @export

constrainLegacy = function(capture.tree = NULL,
                           legacy.taxa = NULL,
                           output.file = NULL,
                           name.match = c("fuzzy", "species", "exact"),
                           drop.unmatched = FALSE,
                           quiet = FALSE) {

  # Debug
  # capture.tree = "capture_iqtree.treefile"
  # legacy.taxa  = names(read_alignment(...))
  # output.file  = "constraint.tre"
  # name.match   = "fuzzy"

  ##############################################################################
  # 0. Parameter checks
  ##############################################################################
  if (is.null(capture.tree)) { stop("Error: capture.tree is required.") }
  if (is.null(legacy.taxa))  { stop("Error: legacy.taxa is required.") }

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
                 error = function(e2) stop("Error: could not parse capture.tree as Newick or NEXUS."))
      }
    )
  } else {
    stop("Error: capture.tree must be a phylo object or a file path string.")
  }

  cap.tips = tree$tip.label
  if (length(cap.tips) == 0) { stop("Error: capture tree has no tip labels.") }

  ##############################################################################
  # 2. Resolve legacy.taxa to a character vector
  ##############################################################################
  if (inherits(legacy.taxa, "phylo")) {
    leg.names = legacy.taxa$tip.label
  } else if (is.character(legacy.taxa)) {
    leg.names = legacy.taxa
  } else {
    stop("Error: legacy.taxa must be a character vector or a phylo object.")
  }

  if (length(leg.names) == 0) { stop("Error: legacy.taxa contains no names.") }

  ##############################################################################
  # 3. Build normalised key function
  ##############################################################################
  if (name.match == "species") {
    # Strip last underscore-delimited token (voucher / specimen ID)
    key.fn = function(x) gsub("_[^_]+$", "", x)
  } else if (name.match == "fuzzy") {
    # Remove all punctuation/separators and lowercase
    key.fn = function(x) tolower(gsub("[-_. ]", "", x))
  } else {
    key.fn = identity
  }

  cap.keys = key.fn(cap.tips)
  leg.keys = key.fn(leg.names)

  ##############################################################################
  # 4. Match each capture tip to exactly one legacy name
  ##############################################################################
  matched.legacy.names = character(length(cap.tips))
  unmatched.tips       = character(0)
  ambiguous.tips       = character(0)

  for (i in seq_along(cap.tips)) {
    hit = which(leg.keys == cap.keys[i])

    if (length(hit) == 1) {
      matched.legacy.names[i] = leg.names[hit]

    } else if (length(hit) == 0) {
      matched.legacy.names[i] = NA_character_
      unmatched.tips = c(unmatched.tips, cap.tips[i])

    } else {
      # Multiple legacy names share the same key.
      # Prefer an exact string match if one exists; otherwise flag as ambiguous.
      exact.hit = which(leg.names[hit] == cap.tips[i])
      if (length(exact.hit) == 1) {
        matched.legacy.names[i] = leg.names[hit[exact.hit]]
      } else {
        matched.legacy.names[i] = leg.names[hit[1]]
        ambiguous.tips = c(ambiguous.tips,
                           paste0("  '", cap.tips[i], "' → ",
                                  paste0("'", leg.names[hit], "'", collapse = ", ")))
      }
    }
  }#end for i

  ##############################################################################
  # 5. Report / handle ambiguous matches
  ##############################################################################
  if (length(ambiguous.tips) > 0 && !quiet) {
    warning(paste0(length(ambiguous.tips),
                   " capture tip(s) matched multiple legacy name(s); ",
                   "the first match was used:\n",
                   paste(ambiguous.tips, collapse = "\n")))
  }

  ##############################################################################
  # 6. Handle unmatched tips
  ##############################################################################
  if (length(unmatched.tips) > 0) {
    unmatched.msg = paste0(
      length(unmatched.tips), " capture tip(s) could not be matched to any ",
      "legacy taxon (name.match = '", name.match, "'):\n",
      paste0("  ", unmatched.tips, collapse = "\n"), "\n",
      "Check for naming discrepancies between the capture tree and your ",
      "integrated alignment.  Set drop.unmatched = TRUE to prune these tips ",
      "from the constraint tree instead of stopping."
    )

    if (drop.unmatched) {
      warning(paste0("Dropping ", length(unmatched.tips),
                     " unmatched tip(s) from constraint tree:\n",
                     paste0("  ", unmatched.tips, collapse = "\n")))
      tree = ape::drop.tip(tree, unmatched.tips)
      keep.idx             = !cap.tips %in% unmatched.tips
      matched.legacy.names = matched.legacy.names[keep.idx]
      cap.tips             = cap.tips[keep.idx]
    } else {
      stop(paste0("Error: ", unmatched.msg))
    }
  }

  ##############################################################################
  # 7. Check for duplicate legacy names after mapping
  #    (two capture tips mapped to the same legacy name would create an invalid tree)
  ##############################################################################
  dup.mapped = duplicated(matched.legacy.names)
  if (any(dup.mapped)) {
    dup.nms = unique(matched.legacy.names[dup.mapped])
    dup.detail = sapply(dup.nms, function(nm) {
      src = cap.tips[matched.legacy.names == nm]
      paste0("  '", nm, "' ← ", paste0("'", src, "'", collapse = ", "))
    })
    stop(paste0(
      "Error: ", length(dup.nms), " legacy name(s) were matched by more than ",
      "one capture tip.  Each legacy taxon can appear only once in a tree:\n",
      paste(dup.detail, collapse = "\n"), "\n",
      "Use a stricter name.match strategy or rename the conflicting taxa."
    ))
  }

  ##############################################################################
  # 8. Rename tips and write output
  ##############################################################################
  tree$tip.label = matched.legacy.names

  if (!is.null(output.file)) {
    ape::write.tree(tree, file = output.file)
    if (!quiet) {
      print(paste0("Constraint tree written to: ", output.file))
    }
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
