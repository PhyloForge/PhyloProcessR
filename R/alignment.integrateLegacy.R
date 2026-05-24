#' @title integrateLegacy
#'
#' @description Integrates legacy (e.g. Sanger or GenBank) sequence alignments into a set
#' of sequence-capture alignments. For each legacy alignment, a consensus sequence is
#' generated and BLASTed against a reference target file to identify the corresponding
#' capture locus. The legacy sequences are then added to the matching capture alignment
#' using MAFFT. Optionally, legacy alignments for loci absent from the capture dataset can
#' be included as stand-alone alignments. When \code{combine.same.sample} is TRUE, samples
#' appearing in both the legacy and capture alignments are merged column-by-column,
#' preferring non-gap characters. Results are written to two output directories:
#' \code{output.directory-only} (legacy-integrated alignments only) and, if
#' \code{include.all.together} is TRUE, \code{output.directory-all} (full merged dataset).
#'
#' @param alignment.directory path to the directory containing the existing sequence-capture
#' alignment files.
#'
#' @param alignment.format format of the sequence-capture alignment files. Accepted values:
#' "phylip" or "fasta".
#'
#' @param output.directory base path for output directories. Two directories are created:
#' \code{output.directory-only} and \code{output.directory-all}.
#'
#' @param legacy.directory path to the directory containing the legacy alignment files to
#' be integrated.
#'
#' @param legacy.format format of the legacy alignment files. Accepted values: "phylip" or
#' "fasta".
#'
#' @param target.markers path to the FASTA file of reference target sequences used to
#' match each legacy alignment to the correct capture locus via BLAST.
#'
#' @param combine.same.sample logical. If TRUE, sequences from the same sample present in
#' both the legacy and capture alignments are merged into a single sequence, preferring
#' non-gap and non-N characters at each site. Default TRUE.
#'
#' @param name.match character. Controls how samples are matched between the legacy and
#' capture alignments when \code{combine.same.sample} is TRUE. Three options:
#' \describe{
#'   \item{\code{"exact"}}{Names must be identical to merge. All legacy sequences are added
#'   to the alignment; those whose names exactly match a capture sequence are merged
#'   column-by-column.}
#'   \item{\code{"species"}}{Strips the trailing specimen/voucher ID (the last
#'   underscore-delimited field) before matching, so that e.g.
#'   \code{Centrolene_bacatum_MZUTI-2436} and \code{Centrolene_bacatum_KU12345} are
#'   treated as the same taxon and merged into a single sequence named
#'   \code{Centrolene_bacatum}. When multiple legacy sequences exist for the same species,
#'   one is pre-selected: the sequence whose full name matches a capture specimen is
#'   preferred; otherwise the most informative (fewest gaps/Ns) sequence is used.}
#'   \item{\code{"fuzzy"}}{Strips all separator characters (hyphens, underscores, dots,
#'   and spaces) and lowercases names before matching, then merges sequences whose
#'   normalised names are identical. This handles common formatting differences including
#'   \code{MZUTI-2436} vs \code{MZUTI_2436} vs \code{MZUTI2436} (separator present,
#'   different separator, or no separator). All legacy sequences that do not match any
#'   capture sequence are added to the alignment as separate rows with their original
#'   names preserved. The merged sequence retains the original capture alignment name.}
#' }
#' Default \code{"exact"}.
#'
#' @param include.uncaptured.legacy logical. If TRUE, legacy alignments for loci not found
#' in the capture dataset are saved to the output as stand-alone alignments. Default FALSE.
#'
#' @param include.all.together logical. If TRUE, all capture alignments are copied to
#' \code{output.directory-all} and updated with legacy-integrated versions where available.
#' Default FALSE.
#'
#' @param threads number of threads passed to BLAST. Default 1.
#'
#' @param memory not currently used; reserved for future parallelisation. Default 1.
#'
#' @param overwrite logical. If TRUE, the output directories are deleted and recreated;
#' if FALSE and they already exist, the function stops with an error. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses BLAST screen output. Default FALSE.
#'
#' @param mafft.path path to the directory containing the MAFFT executable. If NULL, MAFFT
#' is expected to be on the system PATH.
#'
#' @param include.mitochondrial logical. If TRUE, mitochondrial loci absent from the
#' nuclear target marker file can still be integrated by building a second BLAST database
#' from consensus sequences derived from a set of existing mitochondrial capture alignments.
#' When a legacy locus produces no BLAST hit against the nuclear targets, it is re-queried
#' against the mitochondrial database, and if a match is found the corresponding
#' mitochondrial capture alignment is used for integration. Recommended: generate the
#' mitochondrial capture alignments with \href{https://github.com/chutter/MitoTrawlR}{MitoTrawlR}.
#' Default FALSE.
#'
#' @param mito.alignment.directory path to the directory containing the mitochondrial
#' capture alignment files. Only used when \code{include.mitochondrial = TRUE}. The
#' alignment files should cover all mitochondrial loci to be matched.
#'
#' @param mito.alignment.format format of the mitochondrial alignment files. Accepted
#' values: "phylip" or "fasta". Default "phylip".
#'
#' @param blast.path path to the directory containing BLAST executables. If NULL, BLAST
#' tools are expected to be on the system PATH.
#'
#' @return Writes integrated alignment files to \code{output.directory-only} and
#' optionally \code{output.directory-all}. No value is returned to R.
#'
#' @export

integrateLegacy = function(alignment.directory = NULL,
                           alignment.format = "phylip",
                           output.directory = NULL,
                           legacy.directory = NULL,
                           legacy.format = "phylip",
                           target.markers = NULL,
                           combine.same.sample = TRUE,
                           name.match = c("exact", "species", "fuzzy"),
                           include.uncaptured.legacy = FALSE,
                           include.all.together = FALSE,
                           include.mitochondrial = FALSE,
                           mito.alignment.directory = NULL,
                           mito.alignment.format = "phylip",
                           threads = 1,
                           memory = 1,
                           overwrite = FALSE,
                           quiet = FALSE,
                           mafft.path = NULL,
                           blast.path = NULL) {

  # alignment.directory = "/Volumes/LaCie/mitocap_2/Alignments/untrimmed-alignments"
  # alignment.format = "phylip"
  # output.directory = "/Volumes/LaCie/mitocap_2/untrimmed_mt_legacy"
  # legacy.directory = "/Volumes/LaCie/mitocap_2/genbank-legacy"
  # legacy.format = "phylip"
  # target.markers = "/Volumes/LaCie/mitocap_2/reference/refMarkers.fa"
  # combine.same.sample = FALSE
  # include.uncaptured.legacy = FALSE
  # include.all.together = TRUE
  # threads = 10
  # memory = 50
  # overwrite = TRUE
  # quiet = FALSE
  # mafft.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"
  # blast.path = "/Users/chutter/miniconda3/envs/PhyloProcessR/bin"

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  #Same adds to bbmap path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  name.match = match.arg(name.match)

  # Validate mitochondrial inputs
  if (include.mitochondrial == TRUE) {
    if (is.null(mito.alignment.directory)) {
      stop("mito.alignment.directory must be provided when include.mitochondrial = TRUE.")
    }
    if (!dir.exists(mito.alignment.directory)) {
      stop(paste0("mito.alignment.directory not found: ", mito.alignment.directory))
    }
  }

  #Checks this
  if (alignment.directory == output.directory){ stop("You should not overwrite the original alignments.") }

  #Overwrite — check for non-empty directories only; empty dirs are left by interrupted runs
  only.has.files = length(list.files(paste0(output.directory, "-only"))) > 0
  all.has.files  = length(list.files(paste0(output.directory, "-all")))  > 0
  if (only.has.files || all.has.files) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory, "-only"))
      system(paste0("rm -r ", output.directory, "-all"))
    } else {
      stop("Overwrite = FALSE and output directory exists. Either change to TRUE or overwrite manually.")
    }
  }
  dir.create(paste0(output.directory, "-only"), recursive = TRUE, showWarnings = FALSE)
  if (include.all.together == TRUE) { dir.create(paste0(output.directory, "-all"), recursive = TRUE, showWarnings = FALSE) }

  #Gathers alignments
  align.files = list.files(alignment.directory)
  legacy.files = list.files(legacy.directory)

  if (length(align.files) == 0) { stop("alignment files could not be found.") }
  if (length(legacy.files) == 0) { stop("alignment legacy files could not be found.") }

  #headers for the blast db
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  #Make blast database for the nuclear probe loci
  system(paste0(blast.path, "makeblastdb -in ", target.markers,
                " -parse_seqids -dbtype nucl -out target_nucl-blast_db"), ignore.stdout = quiet)

  # Helper: return a sanitized copy of an alignment file if it contains '?' (a valid
  # NEXUS/Sanger missing-data character that Biostrings does not accept).  '?' is
  # replaced with 'N' in the file text before Biostrings reads it.  Returns the
  # original path unchanged if no '?' are found (no temp file created).
  sanitize.align.file = function(file.path) {
    raw.lines = readLines(file.path, warn = FALSE)
    if (!any(grepl("\\?", raw.lines))) { return(list(path = file.path, tmp = FALSE)) }
    san.file = tempfile(fileext = ".tmp.phy")
    writeLines(gsub("\\?", "N", raw.lines), san.file)
    return(list(path = san.file, tmp = TRUE))
  }

  ##########################################################################
  # Build mitochondrial BLAST DB from consensus sequences of mito alignments
  ##########################################################################
  mito.align.files = character(0)
  if (include.mitochondrial == TRUE) {
    print("Building mitochondrial BLAST database from consensus sequences...")
    mito.align.files = list.files(mito.alignment.directory)
    if (length(mito.align.files) == 0) {
      stop("No mitochondrial alignment files found in mito.alignment.directory.")
    }

    mito.cons.list = Biostrings::DNAStringSet()
    for (mf in mito.align.files) {
      mito.locus.name = gsub("\\..*$", "", mf)
      mito.san = sanitize.align.file(paste0(mito.alignment.directory, "/", mf))
      if (mito.alignment.format == "phylip") {
        mito.aln = Biostrings::readDNAMultipleAlignment(file = mito.san$path, format = "phylip")
        mito.aln = Biostrings::DNAStringSet(mito.aln)
      } else {
        mito.aln = Biostrings::readDNAStringSet(mito.san$path)
      }
      if (mito.san$tmp) { file.remove(mito.san$path) }
      # Use most complete sequence as the locus representative in the BLAST DB
      mito.informative = sapply(as.character(mito.aln), function(s) {
        nchar(gsub("[-?nN]", "", s, ignore.case = TRUE))
      })
      mito.con = mito.aln[which.max(mito.informative)]
      names(mito.con) = mito.locus.name
      mito.cons.list = append(mito.cons.list, mito.con)
      rm(mito.aln, mito.con)
    }#end mito files loop

    Biostrings::writeXStringSet(mito.cons.list, filepath = "mito_consensus_references.fa")
    system(paste0(blast.path, "makeblastdb -in mito_consensus_references.fa",
                  " -dbtype nucl -out mito_nucl-blast_db"), ignore.stdout = quiet)
    rm(mito.cons.list)
    print(paste0("Mitochondrial BLAST database built from ", length(mito.align.files), " loci."))
  }#end include.mitochondrial

  ### Reads in the additional stuff
  #add.taxa = Biostrings::readDNAStringSet(sample.markers)
  bait.loci = Biostrings::readDNAStringSet(target.markers)  # loads up fasta file

  # Vectors for summary log: track which output loci carried legacy / mito data
  legacy.loci = character(0)
  mito.loci   = character(0)

  for (i in 1:length(legacy.files)) {

    use.mito = FALSE

    ##############
    #STEP 0: Load in legacy alignment
    ##############
    #Load in alignments
    san = sanitize.align.file(paste0(legacy.directory, "/", legacy.files[i]))
    if (legacy.format == "phylip"){
      align = Biostrings::readDNAMultipleAlignment(file = san$path, format = "phylip")
      align = Biostrings::DNAStringSet(align)
      save.name = gsub("\\..*$", "", legacy.files[i])
    }#end phylip

    if (legacy.format == "fasta"){
      align = Biostrings::readDNAStringSet(san$path)
      save.name = gsub("\\..*$", "", legacy.files[i])
    }#end fasta
    if (san$tmp) { file.remove(san$path) }

    ##############
    #STEP 1: Blast to targets
    ##############

    # Use the most informative (fewest gaps/Ns) single sequence as the BLAST query.
    # A majority-rule consensus of a diverse alignment (e.g. 300 taxa across a whole
    # family) is typically so degenerate it fails to BLAST reliably. A real sequence
    # from the alignment is a much better representative.
    n.informative = sapply(as.character(align), function(s) {
      nchar(gsub("[-?nN]", "", s, ignore.case = TRUE))
    })
    blast.seq = align[which.max(n.informative)]

    #Writes representative sequence to temp file for BLAST
    Biostrings::writeXStringSet(blast.seq, filepath = paste0(save.name, "_query.fa"))

    #Matches samples to loci — use standard blastn for full sensitivity across divergent taxa
    system(paste0(blast.path, "blastn -task blastn -db target_nucl-blast_db -evalue 0.001",
                  " -query ", save.name, "_query.fa -out ", save.name, "_target-blast-match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads), ignore.stdout = quiet, ignore.stderr = quiet)

    #Loads in match data
    match.data = data.table::fread(paste0(save.name, "_target-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    if (nrow(match.data) > 0) { data.table::setnames(match.data, headers) }

    if (nrow(match.data) == 0) {
      # Nuclear BLAST failed — try mitochondrial DB if enabled
      if (include.mitochondrial == TRUE) {
        system(paste0(blast.path, "blastn -task blastn -db mito_nucl-blast_db -evalue 0.001",
                      " -query ", save.name, "_query.fa -out ", save.name, "_mito-blast-match.txt",
                      " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                      " -num_threads ", threads), ignore.stdout = quiet, ignore.stderr = quiet)

        match.data = data.table::fread(paste0(save.name, "_mito-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
        if (nrow(match.data) > 0) { data.table::setnames(match.data, headers) }

        if (nrow(match.data) == 0) {
          print(paste0(save.name, " had no BLAST matches to nuclear or mitochondrial markers. Skipping."))
          system(paste0("rm ", save.name, "*"))
          next
        }
        use.mito = TRUE
      } else {
        print(paste0(save.name, " had no BLAST matches to target markers. Skipping."))
        system(paste0("rm ", save.name, "*"))
        next
      }
    }

    if (nrow(match.data) >= 2){
      match.data = match.data[match.data$bitscore == max(match.data$bitscore),]
    }

    if (nrow(match.data) >= 2){
      stop(paste0(save.name, " matched multiple targets after bitscore filter. Check target file for duplicates."))
    }

    # Search the correct file list depending on whether the match came from mito or nuclear DB
    if (use.mito == TRUE) {
      found.align = mito.align.files[grep(match.data$tName, mito.align.files)]
    } else {
      found.align = align.files[grep(match.data$tName, align.files)]
    }

    if (include.uncaptured.legacy == FALSE){
      if (length(found.align) == 0){
        print(paste0(save.name, " not found in the alignments. Moving to next."))
        next
      }
    }#end if

    if (include.uncaptured.legacy == TRUE){
      if (length(found.align) == 0){
        #Saves them
        write.temp = strsplit(as.character(align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

        #readies for saving
        PhyloProcessR::writePhylip(alignment = aligned.set,
                                   file = paste0(output.directory, "-only/", save.name, ".phy"),
                                   interleave = F,
                                   strict = F)

        print(paste0("No sequence capture alignment found. Included uncaptured ", save.name, " successfully!"))
        next
      }
    }#end if

    ##############
    #STEP 2: Runs MAFFT to add
    ##############
    # Select the source directory and format: mitochondrial or nuclear capture alignments
    src.dir    = if (use.mito == TRUE) mito.alignment.directory else alignment.directory
    src.format = if (use.mito == TRUE) mito.alignment.format    else alignment.format

    #Load in alignments
    old.san = sanitize.align.file(paste0(src.dir, "/", found.align))
    if (src.format == "phylip"){
      old.align = Biostrings::readDNAMultipleAlignment(file = old.san$path, format = "phylip")
      old.align = Biostrings::DNAStringSet(old.align)
      found.name = gsub("\\..*$", "", found.align)
    }#end phylip

    if (src.format == "fasta"){
      old.align = Biostrings::readDNAStringSet(old.san$path)
      found.name = gsub("\\..*$", "", found.align)
    }#end fasta
    if (old.san$tmp) { file.remove(old.san$path) }

    # When matching by species, reduce the legacy alignment to one sequence per species
    # before passing to MAFFT to avoid creating multiple duplicates that complicate merging.
    # Preference: a legacy sequence whose full name matches a capture specimen (museum ID match).
    # Fallback: the sequence with the most informative (non-gap, non-missing) bases.
    if (name.match == "species") {
      stripped.legacy = gsub("_[^_]+$", "", names(align))
      dup.species = unique(stripped.legacy[duplicated(stripped.legacy)])

      if (length(dup.species) > 0) {
        nondup.align  = align[!stripped.legacy %in% dup.species]
        selected.seqs = Biostrings::DNAStringSet()

        for (sp in dup.species) {
          candidates = align[stripped.legacy == sp]

          # Prefer candidate whose full name (with voucher) already exists in capture alignment
          cap.match = candidates[names(candidates) %in% names(old.align)]

          if (length(cap.match) > 0) {
            selected.seqs = append(selected.seqs, cap.match[1])
          } else {
            # No museum ID match — pick the sequence with the most informative bases
            n.informative = sapply(as.character(candidates), function(s) {
              nchar(gsub("[-?nN]", "", s, ignore.case = TRUE))
            })
            selected.seqs = append(selected.seqs, candidates[which.max(n.informative)])
          }
        }#end dup.species loop

        align = append(nondup.align, selected.seqs)
      }
    }# end name.match species pre-selection

    #Aligns and then reverses back to correction orientation
    combo.align = runMafft(sequence.data = old.align,
                           add.contigs = align,
                           save.name = paste0(output.directory, "-only/", found.name),
                           algorithm = "add",
                           adjust.direction = TRUE,
                           threads = 1,
                           cleanup.files = T,
                           quiet = quiet,
                           mafft.path = mafft.path)

    #Checks for failed mafft run
    if (length(combo.align) == 0){ next }
    #Aligns and then reverses back to correction orientation
    names(combo.align) = gsub(pattern = "^_R_", replacement = "", x = names(combo.align))

    # Duplication changes and such
    if (combine.same.sample == TRUE){

      # Identify which sequences in combo.align came from the capture vs legacy source.
      # MAFFT --add outputs newly-added sequences (legacy) first, then the existing
      # alignment (capture) — use position, not names, because when the same specimen
      # appears in both datasets the names are identical and name-matching wrongly
      # classifies the legacy copy as capture, preventing the merge.
      cap.src.names = names(old.align)
      n.legacy      = length(align)
      n.capture     = length(old.align)
      if (length(combo.align) == n.legacy + n.capture) {
        is.cap.seq = c(rep(FALSE, n.legacy), rep(TRUE, n.capture))
      } else {
        # Unexpected output length (MAFFT edge case) — fall back to name matching
        is.cap.seq = names(combo.align) %in% cap.src.names
      }

      # Build normalised keys for matching
      if (name.match == "species") {
        key.fn = function(x) gsub("_[^_]+$", "", x)
      } else if (name.match == "fuzzy") {
        # Strip all separators (hyphens, underscores, dots, spaces) and lowercase,
        # so MZUTI-2436, MZUTI_2436, and MZUTI2436 all normalise to the same key
        key.fn = function(x) tolower(gsub("[-_. ]", "", x))
      } else {
        key.fn = identity
      }

      combo.keys = key.fn(names(combo.align))
      cap.keys   = key.fn(cap.src.names)
      leg.keys   = key.fn(names(align))

      # Only merge keys that appear in BOTH the capture and legacy source sets.
      # This prevents within-dataset duplicates (e.g. two legacy specimens of the
      # same species) from being silently collapsed into one merged sequence.
      dup.keys = intersect(unique(cap.keys), unique(leg.keys))

      if (length(dup.keys) > 0){
        # Start nodup.align with all sequences whose key is not in any cross-match.
        # Unmatched capture-only AND unmatched legacy-only sequences both go here.
        nodup.align = combo.align[!combo.keys %in% dup.keys]

        save.seqs = Biostrings::DNAStringSet()
        for (j in seq_along(dup.keys)){
          dk = dup.keys[j]

          # Separate capture and legacy sequences with this key in combo.align
          cap.idx = which(combo.keys == dk &  is.cap.seq)
          leg.idx = which(combo.keys == dk & !is.cap.seq)

          # Guard: if MAFFT slightly modified a name so we can't find one side, pass through
          if (length(cap.idx) == 0 || length(leg.idx) == 0) {
            nodup.align = append(nodup.align, combo.align[c(cap.idx, leg.idx)])
            next
          }

          # If multiple legacy sequences share this key, keep all but the most informative
          # as separate rows so they are not silently dropped from the output
          if (length(leg.idx) > 1) {
            n.inf = sapply(as.character(combo.align[leg.idx]), function(s)
              nchar(gsub("[-nN?]", "", s, ignore.case = TRUE)))
            best.l  = leg.idx[which.max(n.inf)]
            extras  = leg.idx[leg.idx != best.l]
            nodup.align = append(nodup.align, combo.align[extras])
            leg.idx = best.l
          }

          # Similarly, keep extra capture sequences as-is if more than one match
          if (length(cap.idx) > 1) {
            nodup.align = append(nodup.align, combo.align[cap.idx[-1]])
            cap.idx = cap.idx[1]
          }

          # Merge the one capture + one best-matching legacy sequence column-by-column
          dup.sample = combo.align[c(cap.idx, leg.idx)]
          new.align  = strsplit(as.character(dup.sample), "")
          align.in   = matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = TRUE)

          save.string = character(ncol(align.in))
          for (k in seq_len(ncol(align.in))){
            temp.col = align.in[, k]
            if (length(unique(temp.col)) == 1){
              save.string[k] = temp.col[1]
            } else {
              temp.col = temp.col[!temp.col %in% c("-", "N", "?")]
              save.string[k] = if (length(temp.col) == 0) "-" else temp.col[1]
            }
          }#end k loop

          out.consensus = Biostrings::DNAStringSet(paste0(save.string, collapse = ""))
          # "species" mode: use the collapsed species key; other modes: keep capture name
          names(out.consensus) = if (name.match == "species") dk else names(combo.align)[cap.idx]
          save.seqs = append(save.seqs, out.consensus)

        }#end j

        combo.align = append(save.seqs, nodup.align)
      }#end if dup.keys > 0

    }#end combine.same.sample if

    #Saves them
    write.temp = strsplit(as.character(combo.align), "")
    aligned.set = as.matrix(ape::as.DNAbin(write.temp) )

    #readies for saving
    # Use found.name (the capture locus name) so the file in -only/ correctly
    # overwrites the matching capture alignment when copied to -all/.
    # Using save.name (the legacy locus name) here would leave the original
    # capture-only alignment in -all/ under a different filename, causing
    # concatenateGenes to build gene alignments without legacy samples.
    PhyloProcessR::writePhylip(alignment = aligned.set,
                               file = paste0(output.directory, "-only/", found.name, ".phy"),
                               interleave = F,
                               strict = F)

    print(paste0("Finished ", save.name, " → ", found.name, " legacy integration successfully!"))
    legacy.loci = c(legacy.loci, found.name)
    if (use.mito == TRUE) { mito.loci = c(mito.loci, found.name) }
    system(paste0("rm ", save.name, "*"))
    rm(align, old.align, combo.align, aligned.set, blast.seq, match.data)
    gc()


  }#end loop

  system(paste0("rm target_nucl-blast_db*"))
  if (include.mitochondrial == TRUE) {
    system(paste0("rm mito_nucl-blast_db* mito_consensus_references.fa"))
  }

  ####################################################################################
  if (include.all.together == TRUE){
    # Copy all original capture alignments first.
    # When name.match = "species", rename sequences to collapsed species names
    # so they are consistent with the integrated alignments.

    # Helper: copy/rename a single alignment to the -all directory
    copy.align.to.all = function(src.file, src.directory, src.fmt) {
      out.name = paste0(output.directory, "-all/", gsub("\\..*$", "", src.file), ".phy")
      if (name.match == "species") {
        if (src.fmt == "phylip") {
          cap.align = Biostrings::readDNAMultipleAlignment(
            file = paste0(src.directory, "/", src.file), format = "phylip")
          cap.align = Biostrings::DNAStringSet(cap.align)
        } else {
          cap.align = Biostrings::readDNAStringSet(paste0(src.directory, "/", src.file))
        }
        names(cap.align) = gsub("_[^_]+$", "", names(cap.align))
        write.temp = strsplit(as.character(cap.align), "")
        aligned.set = as.matrix(ape::as.DNAbin(write.temp))
        PhyloProcessR::writePhylip(alignment = aligned.set, file = out.name,
                                   interleave = F, strict = F)
      } else {
        system(paste0("cp ", src.directory, "/", src.file, " ", out.name))
      }
    }#end helper

    # Copy nuclear capture alignments
    for (i in 1:length(align.files)){
      copy.align.to.all(align.files[i], alignment.directory, alignment.format)
    }#end i

    # Copy mitochondrial capture alignments (if enabled)
    if (include.mitochondrial == TRUE && length(mito.align.files) > 0) {
      for (i in 1:length(mito.align.files)){
        copy.align.to.all(mito.align.files[i], mito.alignment.directory, mito.alignment.format)
      }#end i
    }

    # Overwrite with legacy-integrated versions where available (cp overwrites existing)
    new.files = list.files(paste0(output.directory, "-only"))
    for (i in 1:length(new.files)){
      system(paste0("cp ", output.directory, "-only/", new.files[i], " ",
                    output.directory, "-all/", new.files[i]))
    }#end i loop
  }#end if

  ####################################################################################
  ####################################################################################
  ## Summary log
  ####################################################################################
  print("Writing integration summary log...")

  # Choose which directory to summarise:
  # -all includes capture-only loci and is the most complete view;
  # -only is used when include.all.together = FALSE.
  if (include.all.together == TRUE) {
    scan.dir = paste0(output.directory, "-all")
  } else {
    scan.dir = paste0(output.directory, "-only")
  }

  scan.files = list.files(scan.dir)
  only.names = gsub("\\..*$", "", list.files(paste0(output.directory, "-only")))

  # Build per-taxon counts by reading each output alignment
  taxa.counts = list()

  for (sf in scan.files) {
    locus.name = gsub("\\..*$", "", sf)
    aln = Biostrings::readDNAMultipleAlignment(
      file   = paste0(scan.dir, "/", sf),
      format = "phylip")
    aln      = Biostrings::DNAStringSet(aln)
    aln.char = as.character(aln)

    is.mito   = locus.name %in% mito.loci
    is.legacy = locus.name %in% only.names & !is.mito

    for (tx in names(aln)) {
      if (is.null(taxa.counts[[tx]])) {
        taxa.counts[[tx]] = c(Total_Loci          = 0L,
                              Capture_Loci        = 0L,
                              Legacy_Nuclear_Loci = 0L,
                              Legacy_Mito_Loci    = 0L,
                              Total_BP            = 0L)
      }
      bp = nchar(gsub("[-nN?]", "", aln.char[[tx]], ignore.case = TRUE))
      taxa.counts[[tx]]["Total_Loci"] = taxa.counts[[tx]]["Total_Loci"] + 1L
      taxa.counts[[tx]]["Total_BP"]   = taxa.counts[[tx]]["Total_BP"]   + as.integer(bp)

      if (is.mito) {
        taxa.counts[[tx]]["Legacy_Mito_Loci"]    = taxa.counts[[tx]]["Legacy_Mito_Loci"]    + 1L
      } else if (is.legacy) {
        taxa.counts[[tx]]["Legacy_Nuclear_Loci"] = taxa.counts[[tx]]["Legacy_Nuclear_Loci"] + 1L
      } else {
        taxa.counts[[tx]]["Capture_Loci"]        = taxa.counts[[tx]]["Capture_Loci"]        + 1L
      }
    }#end taxon loop
  }#end file loop

  # Assemble data frame
  summary.df = do.call(rbind, lapply(names(taxa.counts), function(tx) {
    r = taxa.counts[[tx]]
    data.frame(
      Taxon               = tx,
      Total_Loci          = r["Total_Loci"],
      Capture_Loci        = r["Capture_Loci"],
      Legacy_Nuclear_Loci = r["Legacy_Nuclear_Loci"],
      Legacy_Mito_Loci    = r["Legacy_Mito_Loci"],
      Legacy_Total        = r["Legacy_Nuclear_Loci"] + r["Legacy_Mito_Loci"],
      Total_BP            = r["Total_BP"],
      Pct_Loci            = round((r["Total_Loci"] / length(scan.files)) * 100, 1),
      stringsAsFactors    = FALSE,
      row.names           = NULL
    )
  }))

  summary.df = summary.df[order(-summary.df$Total_Loci, summary.df$Taxon), ]
  rownames(summary.df) = NULL

  log.file = paste0(output.directory, "-integration_summary.txt")
  write.table(summary.df, file = log.file, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste0("Integration summary written to: ", log.file))

  ####################################################################################

}#end fuction

#END SCRIPT
