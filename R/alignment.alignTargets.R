#' @title alignTargets
#'
#' @description Aligns per-locus sequences extracted from a combined FASTA file using MAFFT,
#' guided by a reference target sequence for each locus. Sequences that diverge from the
#' reference by more than \code{removal.threshold} (pairwise distance) are removed and the
#' locus is realigned. Loci with too few sequences or with duplicate sample names are
#' skipped. Alignments are saved in phylip format, with the reference sequence excluded
#' from the final output. A fractional subset of loci can be processed for parallelisation
#' across jobs.
#'
#' @param targets.to.align path to a FASTA file containing all sequences to be aligned,
#' with names in the format \code{locusName_|_sampleName}.
#'
#' @param target.file path to a FASTA file of reference target sequences, one per locus,
#' used as a guide during alignment and for divergence filtering.
#'
#' @param additional.sequence.directory optional path to an
#' code{extractGenomeTarget()} output directory. Genome target FASTA files are
#' found recursively and added to code{targets.to.align} before sequences are
#' grouped and aligned by locus. Default NULL.
#'
#' @param output.directory path to the directory where phylip alignments will be saved.
#' Default "alignments".
#'
#' @param algorithm MAFFT alignment algorithm to use. Accepted values: "localpair" or
#' "globalpair". Default "localpair".
#'
#' @param min.taxa threshold for the minimum number of sequences required at a locus.
#' Loci with this many or fewer sequences are skipped. Default 4.
#'
#' @param removal.threshold maximum pairwise distance from the reference sequence allowed
#' for a sample sequence to be retained. Sequences exceeding this threshold are removed
#' before a second alignment pass. Default 0.35.
#'
#' @param subset.start lower fractional boundary (0 to 1) in the sorted complete
#' locus list. The lower boundary is exclusive except at zero. Default 0.
#'
#' @param subset.end upper fractional boundary (0 to 1) in the sorted complete
#' locus list. The upper boundary is inclusive. Default 1 (process all loci).
#'
#' @param adjust.direction logical. If TRUE, MAFFT adjusts sequence direction before
#' aligning. Default TRUE.
#'
#' @param remove.reverse.tag logical. If TRUE, the leading \code{_R_} tag added by MAFFT
#' to reverse-complemented sequences is stripped from sequence names. Default TRUE.
#'
#' @param threads number of threads passed to MAFFT. Default 1.
#'
#' @param memory not currently used; reserved for future parallelisation. Default 1.
#'
#' @param overwrite logical. If TRUE, previously completed alignments are overwritten.
#' A partial subset replaces only that subset's alignment files. If FALSE, valid
#' completed alignments are skipped. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses MAFFT screen output. Default TRUE.
#'
#' @param mafft.path path to the directory containing the MAFFT executable. If NULL, MAFFT
#' is expected to be on the system PATH.
#'
#' @return Writes phylip alignment files and dataset-scoped summary CSV files under
#' \code{output.directory}. Partial subsets receive distinct summary-part files.
#' No value is returned to R.
#'
#' @export

alignTargets = function(targets.to.align = NULL,
                        target.file = NULL,
                        output.directory = "alignments",
                        algorithm = c("localpair", "globalpair"),
                        min.taxa = 4,
                        removal.threshold = 0.35,
                        subset.start = 0,
                        subset.end = 1,
                        adjust.direction = TRUE,
                        remove.reverse.tag = TRUE,
                        threads = 1,
                        memory = 1,
                        overwrite = FALSE,
                        quiet = TRUE,
                        mafft.path = NULL,
                        additional.sequence.directory = NULL) {

  # #Debug setup
  # setwd("/Users/chutter/Dropbox/SharewithCarl")
  # targets.to.align = "Venom-Markers-Nov23_to-align.fa"
  # target.file = "venom_loci_updated_Mar12_cdhit95_duplicate_exons_renamed_Feb2023_FINAL.fa"
  # output.directory = "alignments"
  #
  # #Main settings
  # subset.start = 0
  # subset.end = 1
  # min.taxa = 4
  # threads = 4
  # memory = 8
  # overwrite = TRUE
  # resume = FALSE
  # quiet = TRUE
  # adjust.direction = TRUE
  # algorithm = "localpair"
  # remove.reverse.tag = TRUE
  # removal.threshold = 0.35
  #
  # #program paths
  # mafft.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  algorithm = match.arg(algorithm)

  valid.fraction = function(x) is.numeric(x) && length(x) == 1 && is.finite(x)
  if (!valid.fraction(subset.start) || !valid.fraction(subset.end) ||
      subset.start < 0 || subset.end > 1 || subset.start > subset.end) {
    stop("subset.start and subset.end must be finite scalars with 0 <= start <= end <= 1.")
  }

  #Same adds to bbmap path
  if (is.null(mafft.path) == FALSE){
    b.string = unlist(strsplit(mafft.path, ""))
    if (b.string[length(b.string)] != "/") {
      mafft.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { mafft.path = "" }

  #Initial checks
  if (is.null(targets.to.align) == T){ stop("A fasta file of targets is needed for alignment.") }
  if (!file.exists(targets.to.align)) stop("Targets-to-align FASTA file not found.")
  if (is.null(target.file) || !file.exists(target.file)) stop("Reference target FASTA file not found.")

  input.paths = normalizePath(c(targets.to.align, target.file), mustWork = TRUE)
  output.path = normalizePath(output.directory, mustWork = FALSE)
  if (any(dirname(input.paths) == output.path)) {
    stop("output.directory contains an input FASTA and cannot be safely overwritten.")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE && subset.start == 0 && subset.end == 1){
      unlink(output.directory, recursive = TRUE)
      dir.create(output.directory, recursive = TRUE)
    }
  } else { dir.create(output.directory, recursive = TRUE) }

  #Load sample file from probe matching step
  all.data = Biostrings::readDNAStringSet(targets.to.align)   # loads up fasta file

  if (!is.null(additional.sequence.directory)) {
    if (!dir.exists(additional.sequence.directory)) {
      stop("Additional sequence directory not found: ",
           additional.sequence.directory)
    }

    additional.files = list.files(
      additional.sequence.directory,
      pattern = "_target-matches[.]fa(sta)?$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
    if (length(additional.files) == 0) {
      stop("No extractGenomeTarget FASTA files found in: ",
           additional.sequence.directory)
    }

    additional.data = Biostrings::DNAStringSet()
    for (additional.file in additional.files) {
      additional.data = append(
        additional.data,
        Biostrings::readDNAStringSet(additional.file)
      )
    }

    duplicate.ids = intersect(names(all.data), names(additional.data))
    duplicate.ids = union(
      duplicate.ids,
      unique(names(additional.data)[duplicated(names(additional.data))])
    )
    if (length(duplicate.ids) > 0) {
      stop(
        "Genome target files contain sequence IDs that are already in use: ",
        paste(duplicate.ids, collapse = ", ")
      )
    }
    all.data = append(all.data, additional.data)
  }

  sequence.loci = sub("_\\|_.*$", "", names(all.data))
  locus.names = sort(unique(sequence.loci))

  target.seqs = Biostrings::readDNAStringSet(target.file)   # loads up fasta file
  reference.loci = sub("_\\|_.*$", "", names(target.seqs))

  first = floor(subset.start * length(locus.names)) + 1L
  last = floor(subset.end * length(locus.names))
  if (first > last || length(locus.names) == 0) return(invisible(NULL))
  locus.names = locus.names[seq.int(first, last)]

  if (overwrite && !(subset.start == 0 && subset.end == 1)) {
    unlink(file.path(output.directory, paste0(locus.names, ".phy")))
  }

  #Checks for alignments already done and removes from the to-do list
  if (overwrite == FALSE){
    done = list.files(output.directory, pattern = "\\.phy$", full.names = FALSE)
    valid = vapply(file.path(output.directory, done), function(path) {
      tryCatch({
        alignment = ape::read.dna(path, format = "sequential")
        nrow(as.matrix(alignment)) > 0 && ncol(as.matrix(alignment)) > 0
      }, error = function(e) FALSE)
    }, logical(1))
    if (any(!valid)) {
      print(paste(sum(!valid), "malformed alignment output(s) will be rebuilt."))
      unlink(file.path(output.directory, done[!valid]))
      done = done[valid]
    }
    locus.names = locus.names[!locus.names %in% sub("\\.phy$", "", done)]
  }
  if (length(locus.names) == 0){ return(invisible(NULL)) }

  # Log tracking
  log.directory = file.path(output.directory, "logs")
  dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
  subset.tag = if (subset.start == 0 && subset.end == 1) "" else
    paste0("_part_", format(subset.start, trim = TRUE), "_",
           format(subset.end, trim = TRUE))
  locus.log    = vector("list", length(locus.names))
  aligned.samples = character(0)

  #Loops through each locus and writes each species to end of file
  for (i in seq_along(locus.names)) {

    #Match probe names to contig names to acquire data
    match.data = all.data[sequence.loci == locus.names[i]]

    #STEP 1: Throw out loci if there are too few taxa
    if (length(names(match.data)) < min.taxa){
      print(paste0(locus.names[i], " had too few taxa"))
      locus.log[[i]] = data.frame(
        Locus = locus.names[i], N_taxa_initial = length(match.data),
        N_taxa_aligned = 0L, Alignment_length = NA_integer_,
        Mean_pct_missing = NA_real_, Status = "too_few_taxa",
        stringsAsFactors = FALSE)
      next
    }

    #STEP 2: Sets up fasta for aligning
    names(match.data) = gsub(pattern = ".*_\\|_", replacement = "", x = names(match.data))

    #Gets reference locus
    ref.locus = target.seqs[reference.loci == locus.names[i]]
    if (length(ref.locus) == 0) {
      print(paste0(locus.names[i], ": no reference sequence found in target file -- skipping."))
      locus.log[[i]] = data.frame(
        Locus = locus.names[i], N_taxa_initial = length(match.data),
        N_taxa_aligned = 0L, Alignment_length = NA_integer_,
        Mean_pct_missing = NA_real_, Status = "no_reference",
        stringsAsFactors = FALSE)
      next
    }
    if (length(ref.locus) > 1) {
      print(paste0(locus.names[i], ": multiple reference sequences found -- skipping."))
      locus.log[[i]] = data.frame(
        Locus = locus.names[i], N_taxa_initial = length(match.data),
        N_taxa_aligned = 0L, Alignment_length = NA_integer_,
        Mean_pct_missing = NA_real_, Status = "duplicate_reference",
        stringsAsFactors = FALSE)
      next
    }
    names(ref.locus) = paste("Reference_Locus")
    final.loci = append(match.data, ref.locus)

    # Checks for duplicates
    dup.names = match.data[duplicated(names(match.data)), ]
    if (length(dup.names) != 0) {
      print(paste0(locus.names[i], " did not successfully align. Duplicate samples, likely paralog."))
      locus.log[[i]] = data.frame(
        Locus = locus.names[i], N_taxa_initial = length(match.data),
        N_taxa_aligned = 0L, Alignment_length = NA_integer_,
        Mean_pct_missing = NA_real_, Status = "duplicate_samples",
        stringsAsFactors = FALSE)
      next
    } # end if

    #STEP 3: Runs MAFFT to align
    #Aligns and then reverses back to correction orientation
    alignment = runMafft(sequence.data = final.loci,
                         save.name = paste0(output.directory, "/", locus.names[i]),
                         algorithm = algorithm,
                         adjust.direction = adjust.direction,
                         threads = threads,
                         cleanup.files = T,
                         quiet = quiet,
                         mafft.path = mafft.path)

    #Checks for failed mafft run
    if (length(alignment) == 0){
      print(paste0(locus.names[i], " did not successfully align."))
      locus.log[[i]] = data.frame(
        Locus = locus.names[i], N_taxa_initial = length(match.data),
        N_taxa_aligned = 0L, Alignment_length = NA_integer_,
        Mean_pct_missing = NA_real_, Status = "mafft_failed",
        stringsAsFactors = FALSE)
      next }

    #Removes the reverse name
    if (remove.reverse.tag == TRUE){
      names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
    }#end if


    #Gets the divergence to make sure not crazy
    diff = pairwiseDistanceTarget(alignment, "Reference_Locus")
    undefined.seqs = names(diff)[is.na(diff)]
    if (length(undefined.seqs) > 0) {
      print(paste0(locus.names[i], ": removing ", length(undefined.seqs),
                   " sequence(s) with no comparable reference positions."))
    }
    bad.seqs = names(diff)[which(is.na(diff) | diff >= removal.threshold)]
    rem.align = alignment[!names(alignment) %in% bad.seqs]

    # Moves onto next loop in there are no good sequences
    if (length(rem.align) <= as.numeric(min.taxa)){
      #Deletes old files
      print(paste(locus.names[i], " had too few taxa", sep = ""))
      locus.log[[i]] = data.frame(
        Locus = locus.names[i], N_taxa_initial = length(match.data),
        N_taxa_aligned = 0L, Alignment_length = NA_integer_,
        Mean_pct_missing = NA_real_, Status = "too_few_taxa_after_filter",
        stringsAsFactors = FALSE)
      next }

    ### realign if bad seqs removed
    if (length(bad.seqs) != 0){

      alignment = runMafft(sequence.data = rem.align,
                           save.name = paste0(output.directory, "/", locus.names[i]),
                           algorithm = algorithm,
                           adjust.direction = adjust.direction,
                           threads = threads,
                           cleanup.files = T,
                           quiet = quiet,
                           mafft.path = mafft.path)

      #Checks for failed mafft run
      if (length(alignment) == 0){
        print(paste0(locus.names[i], " did not successfully align."))
        locus.log[[i]] = data.frame(
          Locus = locus.names[i], N_taxa_initial = length(match.data),
          N_taxa_aligned = 0L, Alignment_length = NA_integer_,
          Mean_pct_missing = NA_real_, Status = "mafft_failed_realign",
          stringsAsFactors = FALSE)
        next }

      #Removes the reverse name
      if (remove.reverse.tag == TRUE){
        names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))
      }#end if

    } # end bad.seqs if

    #Saves prelim exon file
    red.align = alignment[!names(alignment) %in% "Reference_Locus"]
    new.align = strsplit(as.character(red.align), "")
    aligned.set = as.matrix(ape::as.DNAbin(new.align))

    # Compute per-locus stats for log
    aln.mat    = as.character(aligned.set)
    aln.len    = ncol(aligned.set)
    gap.counts = rowSums(aln.mat == "-" | aln.mat == "n" | aln.mat == "N")
    pct.miss   = round(mean(gap.counts / aln.len) * 100, 2)

    locus.log[[i]] = data.frame(
      Locus            = locus.names[i],
      N_taxa_initial   = length(match.data),
      N_taxa_aligned   = length(red.align),
      Alignment_length = aln.len,
      Mean_pct_missing = pct.miss,
      Status           = "aligned",
      stringsAsFactors = FALSE)

    aligned.samples = c(aligned.samples, rownames(aligned.set))

    #readies for saving
    final.file = file.path(output.directory, paste0(locus.names[i], ".phy"))
    temp.file = tempfile(paste0(locus.names[i], "-"), tmpdir = output.directory)
    PhyloProcessR::writePhylip(alignment = aligned.set,
                               file = temp.file,
                               interleave = F,
                               strict = F)
    if (!file.rename(temp.file, final.file)) {
      unlink(temp.file)
      stop("Could not publish completed alignment: ", final.file)
    }

    #Deletes old files
    print(paste0(locus.names[i], " alignment saved."))

  }# end big i loop

  ##########################################################################
  # Write alignment logs
  ##########################################################################
  log.df = do.call(rbind, locus.log[!sapply(locus.log, is.null)])

  if (!is.null(log.df) && nrow(log.df) > 0) {

    # Append to any existing log (handles subset / resume runs)
    locus.log.file = file.path(log.directory,
                               paste0("alignTargets_locus_summary", subset.tag, ".csv"))
    if (file.exists(locus.log.file)) {
      existing = read.csv(locus.log.file, stringsAsFactors = FALSE)
      existing = existing[!existing$Locus %in% log.df$Locus, ]
      log.df   = rbind(existing, log.df)
    }
    write.csv(log.df, file = locus.log.file, row.names = FALSE)
    print(paste0("Per-locus alignment summary: ", locus.log.file))
  }

  completed.files = list.files(output.directory, pattern = "\\.phy$", full.names = TRUE)
  completed.samples = unlist(lapply(completed.files, function(path) {
    tryCatch(rownames(as.matrix(ape::read.dna(path, format = "sequential"))),
             error = function(e) character())
  }), use.names = FALSE)
  if (length(completed.samples) > 0) {
    samp.tbl  = as.data.frame(table(completed.samples), stringsAsFactors = FALSE)
    names(samp.tbl) = c("Sample", "N_loci_aligned")
    samp.tbl  = samp.tbl[order(-samp.tbl$N_loci_aligned), ]

    samp.log.file = file.path(log.directory,
                              paste0("alignTargets_sample_summary", subset.tag, ".csv"))
    write.csv(samp.tbl, file = samp.log.file, row.names = FALSE)
    print(paste0("Per-sample alignment summary: ", samp.log.file))
  }

}# end function


#END SCRIPT


#END SCRIPT
