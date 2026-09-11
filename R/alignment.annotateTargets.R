.summarizeTargetCandidate = function(hits = NULL, sequence = NULL) {

  hits = hits[order(-hits$bitscore, -hits$matches, hits$evalue), ]
  target.length = max(hits$qLen)
  covered = rep(FALSE, target.length)
  supported.bp = 0
  identity.total = 0
  score.total = 0

  for (row in seq_len(nrow(hits))) {
    start = max(1, min(hits$qStart[row], hits$qEnd[row]))
    end = min(target.length, max(hits$qStart[row], hits$qEnd[row]))
    positions = seq.int(start, end)
    new.bp = sum(!covered[positions])
    if (new.bp == 0) next

    covered[positions] = TRUE
    supported.bp = supported.bp + new.bp
    identity.total = identity.total + (new.bp * hits$pident[row])
    score.total = score.total +
      (hits$bitscore[row] * new.bp / max(1, hits$matches[row]))
  }

  sequence.text = toupper(as.character(sequence))
  sequence.length = nchar(sequence.text)
  non.n.length = nchar(gsub("N", "", sequence.text, fixed = TRUE))
  iupac.count = nchar(gsub("[^RYKMSWBDHV]", "", sequence.text))

  data.frame(
    Search_hit_count = nrow(hits),
    Supported_bp = supported.bp,
    Target_length = target.length,
    Target_coverage = supported.bp / target.length,
    Weighted_identity = identity.total / supported.bp,
    Total_bitscore = score.total,
    Best_evalue = min(hits$evalue),
    Contig_length = sequence.length,
    Non_N_length = non.n.length,
    IUPAC_count = iupac.count,
    IUPAC_proportion = iupac.count / sequence.length,
    stringsAsFactors = FALSE
  )
}

#' @title annotateTargets
#'
#' @description Annotates assembly contigs. It matches them to a set of target
#' marker sequences. It searches curated contigs against the target file without
#' collapsing similar copies.
#' \code{search.method} selects the program. The default is LAST. The function
#' combines overlapping search hits into one candidate per target and contig,
#' then filters candidates by identity, supported length, and coverage. It ranks
#' candidates and applies \code{paralog.action}. The annotated contigs for each
#' sample are saved as a per-sample FASTA file in \code{output.directory}. A combined
#' FASTA file suitable for downstream alignment (named
#' \code{alignment.contig.name_to-align.fa}) and a summary CSV are written to
#' the working directory.
#'
#' @details The structural curation moved to \code{curateTargetContigs}, which
#'   runs at the end of workflow 2. That function joins the fragments of one
#'   target and cuts apart a contig that spans more than one target, before the
#'   variant caller maps reads to the contigs. This function keeps the paralog
#'   policy, the sample naming, and the alignment output. Give it contigs that
#'   \code{curateTargetContigs} has already curated. A secondary candidate is
#'   competitive only when it passes the absolute filters and the configured
#'   score, coverage, and identity comparisons with the best candidate.
#'
#' @param assembly.directory path to the directory containing per-sample contig FASTA files
#' (one file per sample, named \code{sampleName.fa}).
#'
#' @param target.file path to the FASTA file of target marker sequences used for the search
#' matching.
#'
#' @param alignment.contig.name base name (without extension) used for the combined output
#' FASTA and summary CSV files. Default "annotated-contigs-all".
#'
#' @param output.directory path to the directory where per-sample annotated contig files
#' will be saved. Default "annotated-contigs".
#'
#' @param min.match.percent minimum percent identity required to retain a hit.
#' Default 60.
#'
#' @param min.match.length minimum alignment length (in bp) required to retain a hit.
#' The N padding that joins two fragments of one target is not counted. Default 50.
#'
#' @param min.match.coverage minimum proportion of the target sequence length that must be
#' covered by the hit (expressed as a percentage). Default 30.
#'
#' @param paralog.action how competing target copies are handled. \code{"exclude"}
#' removes all copies when a second candidate is close to the best candidate.
#' \code{"best"} always keeps the best candidate. Default \code{"exclude"}.
#'
#' @param paralog.score.ratio minimum total-score ratio between the second and
#' best candidates for the second candidate to be considered competitive.
#' Default 0.8.
#'
#' @param paralog.coverage.ratio minimum target-coverage ratio between the second
#' and best candidates for the second candidate to be considered competitive.
#' Default 0.8.
#'
#' @param paralog.identity.delta maximum weighted-identity difference, in
#' percentage points, between the best and second candidates. Default 5.
#'
#' @param paralog.directory directory for all qualifying copies from targets
#' with more than one candidate. Default \code{"paralog-contigs"}.
#'
#' @param threads number of parallel threads to use. Default 1.
#'
#' @param memory total memory (in GB) to allocate across all threads. Default 1.
#'
#' @param search.method which program matches the target markers to the contigs.
#'   \code{"last"} (default) uses LAST, which matches a contig that is up to
#'   about 35 percent divergent from its target. \code{"blast"} uses the
#'   previous \code{blastn dc-megablast} search, which loses a contig past about
#'   25 percent divergence. Numbers in HANDOFF.md.
#'
#' @param last.path path to the directory containing \code{lastdb} and
#'   \code{lastal}. Only needed when \code{search.method = "last"}. If NULL the
#'   programs must be on the system PATH.
#'
#' @param blast.path path to the directory containing BLAST executables. Only
#'   needed when \code{search.method = "blast"}. If NULL, BLAST
#' tools are expected to be on the system PATH.
#'
#' @param overwrite logical. If TRUE, previously completed samples are reprocessed;
#' if FALSE, they are skipped. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses the search program screen output. Default TRUE.
#'
#' @return Writes per-sample annotated FASTA files to \code{output.directory}, a combined
#' FASTA file for alignment, qualifying multi-copy contigs to
#' \code{paralog.directory}, and summary CSV files. These logs are also written:
#' \itemize{
#'   \item \code{logs/sample_logs/<Sample>_target-candidates.csv} -- one row per
#'     target-contig candidate with absolute filters, relative evidence, and decision.
#'   \item \code{logs/sample_logs/<Sample>_target-hits-raw.csv} -- the original
#'     search-hit rows used to calculate candidate evidence.
#'   \item \code{logs/annotateTargets_summary.csv} -- one row per sample summarising
#'     input contigs, retained targets, ambiguous targets, and paralog candidates.
#' }
#' No value is returned to R.
#'
#' @export

annotateTargets = function(assembly.directory = NULL,
                            target.file = NULL,
                            alignment.contig.name = "annotated-contigs-all",
                            output.directory = "annotated-contigs",
                            min.match.percent = 60,
                            min.match.length = 50,
                            min.match.coverage = 30,
                            paralog.action = c("exclude", "best"),
                            paralog.score.ratio = 0.8,
                            paralog.coverage.ratio = 0.8,
                            paralog.identity.delta = 5,
                            paralog.directory = "paralog-contigs",
                            threads = 1,
                            memory = 1,
                            blast.path = NULL,
                            last.path = NULL,
                            search.method = c("last", "blast"),
                            overwrite = FALSE,
                            quiet = TRUE) {

  search.method = match.arg(search.method)
  paralog.action = match.arg(paralog.action)
  blast.path = .programPrefix(blast.path)
  last.path = .programPrefix(last.path)

  ratio.settings = c(paralog.score.ratio, paralog.coverage.ratio)
  if (any(!is.finite(ratio.settings)) || any(ratio.settings < 0) ||
      any(ratio.settings > 1)) {
    stop("Paralog score and coverage ratios must be between 0 and 1.")
  }
  if (length(paralog.identity.delta) != 1 ||
      !is.finite(paralog.identity.delta) || paralog.identity.delta < 0) {
    stop("paralog.identity.delta must be one non-negative finite value.")
  }

  #Initial checks
  if (is.null(assembly.directory) || !dir.exists(assembly.directory)) {
    stop("Assembly directory not found. Please check the path.")
  }
  if (is.null(target.file) == TRUE){ stop("A fasta file of targets to match to assembly contigs is needed.") }
  if (file.exists(target.file) == FALSE){ stop("Target file not found. Please check path / use full path.") }
  assembly.path = normalizePath(assembly.directory, mustWork = TRUE)
  output.path = normalizePath(output.directory, mustWork = FALSE)
  paralog.path = normalizePath(paralog.directory, mustWork = FALSE)
  target.directory = dirname(normalizePath(target.file))
  if (assembly.path %in% c(output.path, paralog.path) ||
      output.path == paralog.path ||
      target.directory %in% c(output.path, paralog.path)) {
    stop("You should not overwrite the original contigs or target file.")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      unlink(output.directory, recursive = TRUE)
      dir.create(output.directory, recursive = TRUE)
    }
  } else { dir.create(output.directory, recursive = TRUE) }

  if (dir.exists(paralog.directory) == TRUE) {
    if (overwrite == TRUE) unlink(paralog.directory, recursive = TRUE)
  }
  dir.create(paralog.directory, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists("logs/sample_logs")) {
    dir.create("logs/sample_logs", recursive = TRUE, showWarnings = FALSE)
  }

  #Gets contig file names
  file.names = list.files(assembly.directory, pattern = "\\.fa$", full.names = FALSE)
  if (length(file.names) == 0) stop("No .fa files found in assembly.directory.")

  #headers for the blast db
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  # Process samples in parallel.
  results = parallel::mclapply(seq_along(file.names), function(i) {
  tryCatch({

    #Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])

    #Checks if this has been done already (before creating any directory)
    if (overwrite == FALSE){
      if (file.exists(file.path(output.directory, paste0(sample, ".fa"))) ||
          file.exists(file.path(output.directory, paste0(sample, ".complete-empty")))) {
        print(paste0(sample, " already finished, skipping. Set overwrite = TRUE to redo."))
        return(NULL)
      }
    }#end

    # Temporary working directory kept inside logs so output.directory stays clean
    species.dir = paste0("logs/sample_logs/", sample)
    if (!dir.exists(species.dir)){ dir.create(species.dir, recursive = TRUE, showWarnings = FALSE) }

    # Read the curated contigs without collapsing similar copies.
    input.file = file.path(assembly.directory, file.names[i])
    all.data = Biostrings::readDNAStringSet(input.file, format = "fasta")
    if (length(all.data) == 0) stop(sample, " has an empty input FASTA.")
    source.names = names(all.data)
    names(all.data) = paste0("contig_", seq_along(all.data))

    # Writes the final loci
    final.loci = as.list(as.character(all.data))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(species.dir, "/", sample, "_rename.fa"),
      nbchar = 1000000, as.string = TRUE, open = "w"
    )

    #########################################################################
    #Part B: Blasting
    #########################################################################

    # The database holds the contigs of this sample and the query is the target
    # file, so a hit names the target first and the contig second.
    # One thread per search, because the samples already run in parallel.
    if (search.method == "last") {
      .lastBuildDB(reference.file = paste0(species.dir, "/", sample, "_rename.fa"),
                   db.prefix = paste0(species.dir, "/", sample, "_last_db"),
                   lastdb.command = paste0(last.path, "lastdb"),
                   threads = 1,
                   quiet = quiet)

      .lastSearch(query.file = target.file,
                  db.prefix = paste0(species.dir, "/", sample, "_last_db"),
                  out.file = paste0(species.dir, "/", sample, "_target-blast-match.txt"),
                  lastal.command = paste0(last.path, "lastal"),
                  threads = 1,
                  quiet = quiet)
    } else {
      blast.db = file.path(species.dir, paste0(sample, "_nucl-blast_db"))
      .runCommand(paste0(shQuote(paste0(blast.path, "makeblastdb")), " -in ",
                         shQuote(file.path(species.dir, paste0(sample, "_rename.fa"))),
                         " -parse_seqids -dbtype nucl -out ", shQuote(blast.db)),
                  quiet = quiet, task = paste(sample, "BLAST database"))

      match.file = file.path(species.dir, paste0(sample, "_target-blast-match.txt"))
      unlink(match.file)
      .runCommand(paste0(
        shQuote(paste0(blast.path, "blastn")),
        " -task dc-megablast -db ", shQuote(blast.db), " -evalue 0.001",
        " -query ", shQuote(target.file), " -out ", shQuote(match.file),
        " -outfmt \"6 qseqid sseqid pident length mismatch gapopen ",
        "qstart qend sstart send evalue bitscore qlen slen gaps\" ",
        " -num_threads 1"
      ), quiet = quiet, task = paste(sample, "BLAST search"))
    }

    # Remove the search database and the large intermediate contig files
    unlink(c(list.files(species.dir, pattern = "nucl-blast_db|_last_db", full.names = TRUE),
             file.path(species.dir, paste0(sample, "_rename.fa"))))

    #Loads in match data
    match.file = file.path(species.dir, paste0(sample, "_target-blast-match.txt"))
    if (file.size(match.file) == 0) {
      raw.log = data.frame(matrix(nrow = 0, ncol = length(c("Sample", headers))))
      names(raw.log) = c("Sample", headers)
      write.csv(raw.log,
                file.path("logs/sample_logs", paste0(sample, "_target-hits-raw.csv")),
                row.names = FALSE)
      evidence.columns = c(
        "Sample", "Target", "Contig", "Source_contig", "Candidate_rank",
        "Search_hit_count", "Supported_bp", "Target_length", "Target_coverage",
        "Weighted_identity", "Total_bitscore", "Best_evalue", "Contig_length",
        "Non_N_length", "IUPAC_count", "IUPAC_proportion",
        "Best_score_ratio", "Best_coverage_ratio", "Identity_difference",
        "Passes_absolute_filters", "Competitive_with_best",
        "Selected_for_primary", "Saved_as_paralog", "Decision")
      evidence = data.frame(matrix(nrow = 0, ncol = length(evidence.columns)))
      names(evidence) = evidence.columns
      write.csv(evidence,
                file.path("logs/sample_logs", paste0(sample, "_target-candidates.csv")),
                row.names = FALSE)
      file.create(file.path(output.directory, paste0(sample, ".complete-empty")))
      return(data.frame(
        Sample = sample, InputContigs = length(all.data), TargetsMatched = 0L,
        AnnotatedTargets = 0L, AmbiguousTargets = 0L, ParalogCandidates = 0L,
        MeanPident = NA_real_, MeanBitscore = NA_real_, MaxBitscore = NA_real_,
        stringsAsFactors = FALSE))
    }
    match.data = data.table::fread(match.file, sep = "\t", header = F, stringsAsFactors = FALSE)
    data.table::setnames(match.data, headers)

    raw.log = as.data.frame(match.data)
    raw.log$Sample = sample
    raw.log = raw.log[, c("Sample", headers)]
    write.csv(raw.log,
              file.path("logs/sample_logs", paste0(sample, "_target-hits-raw.csv")),
              row.names = FALSE)

    candidate.keys = unique(match.data[, .(qName, tName)])
    candidate.rows = vector("list", nrow(candidate.keys))
    for (candidate.index in seq_len(nrow(candidate.keys))) {
      target.name = candidate.keys$qName[candidate.index]
      contig.name = candidate.keys$tName[candidate.index]
      hits = match.data[qName == target.name & tName == contig.name]
      candidate = .summarizeTargetCandidate(hits, all.data[contig.name])
      candidate$Sample = sample
      candidate$Target = target.name
      candidate$Contig = contig.name
      candidate$Source_contig = source.names[match(contig.name, names(all.data))]
      candidate.rows[[candidate.index]] = candidate
    }
    candidate.data = data.table::rbindlist(candidate.rows, fill = TRUE)
    candidate.data[, Passes_absolute_filters :=
      Weighted_identity >= min.match.percent &
      Supported_bp >= min.match.length &
      Target_coverage >= (min.match.coverage / 100)]

    candidate.data[, Candidate_rank := NA_integer_]
    candidate.data[, Best_score_ratio := NA_real_]
    candidate.data[, Best_coverage_ratio := NA_real_]
    candidate.data[, Identity_difference := NA_real_]
    candidate.data[, Competitive_with_best := FALSE]
    candidate.data[, Selected_for_primary := FALSE]
    candidate.data[, Saved_as_paralog := FALSE]

    qualifying = candidate.data[Passes_absolute_filters == TRUE]
    ambiguous.targets = character(0)

    for (target.name in unique(qualifying$Target)) {
      rows = which(candidate.data$Target == target.name &
                   candidate.data$Passes_absolute_filters)
      rows = rows[order(-candidate.data$Total_bitscore[rows],
                        -candidate.data$Target_coverage[rows],
                        -candidate.data$Weighted_identity[rows])]

      best.row = rows[1]
      candidate.data$Candidate_rank[rows] = seq_along(rows)
      candidate.data$Best_score_ratio[rows] =
        candidate.data$Total_bitscore[rows] /
        candidate.data$Total_bitscore[best.row]
      candidate.data$Best_coverage_ratio[rows] =
        candidate.data$Target_coverage[rows] /
        candidate.data$Target_coverage[best.row]
      candidate.data$Identity_difference[rows] =
        candidate.data$Weighted_identity[best.row] -
        candidate.data$Weighted_identity[rows]

      if (length(rows) > 1) {
        other.rows = rows[-1]
        competitive =
          candidate.data$Best_score_ratio[other.rows] >= paralog.score.ratio &
          candidate.data$Best_coverage_ratio[other.rows] >= paralog.coverage.ratio &
          candidate.data$Identity_difference[other.rows] <= paralog.identity.delta
        candidate.data$Competitive_with_best[other.rows] = competitive
        candidate.data$Saved_as_paralog[rows] = TRUE
        if (any(competitive)) ambiguous.targets = c(ambiguous.targets, target.name)
      }

      keep.best = paralog.action == "best" ||
        !target.name %in% ambiguous.targets
      if (keep.best) candidate.data$Selected_for_primary[best.row] = TRUE
    }

    candidate.data[, Decision := "weak_match"]
    candidate.data[Passes_absolute_filters == TRUE, Decision := "secondary_candidate"]
    candidate.data[Competitive_with_best == TRUE, Decision := "competitive_paralog"]
    candidate.data[Selected_for_primary == TRUE, Decision := "selected"]
    if (paralog.action == "exclude" && length(ambiguous.targets) > 0) {
      candidate.data[Target %in% ambiguous.targets & Candidate_rank == 1,
                     Decision := "best_excluded_due_to_competitor"]
    }

    evidence.columns = c(
      "Sample", "Target", "Contig", "Source_contig", "Candidate_rank",
      "Search_hit_count", "Supported_bp", "Target_length", "Target_coverage",
      "Weighted_identity", "Total_bitscore", "Best_evalue", "Contig_length",
      "Non_N_length", "IUPAC_count", "IUPAC_proportion",
      "Best_score_ratio", "Best_coverage_ratio", "Identity_difference",
      "Passes_absolute_filters", "Competitive_with_best",
      "Selected_for_primary", "Saved_as_paralog", "Decision")
    evidence = as.data.frame(candidate.data)[, evidence.columns]
    write.csv(evidence,
              file.path("logs/sample_logs", paste0(sample, "_target-candidates.csv")),
              row.names = FALSE)

    selected = candidate.data[Selected_for_primary == TRUE]
    if (nrow(selected) > 0) {
      selected = selected[order(Target)]
      selected.seqs = all.data[selected$Contig]
      names(selected.seqs) = paste0(selected$Target, "_|_", sample)
      final.file = file.path(output.directory, paste0(sample, ".fa"))
      temp.file = tempfile(paste0(sample, "-"), tmpdir = output.directory)
      PhyloProcessR::writeFasta(
        sequences = as.list(as.character(selected.seqs)),
        names = names(selected.seqs), file = temp.file,
        nbchar = 1000000, as.string = TRUE)
      if (!file.rename(temp.file, final.file)) {
        unlink(temp.file)
        stop("Could not publish annotated FASTA for ", sample)
      }
    } else {
      file.create(file.path(output.directory, paste0(sample, ".complete-empty")))
    }

    paralog.candidates = candidate.data[Saved_as_paralog == TRUE]
    if (nrow(paralog.candidates) > 0) {
      paralog.candidates = paralog.candidates[order(Target, Candidate_rank)]
      paralog.seqs = all.data[paralog.candidates$Contig]
      names(paralog.seqs) = paste0(
        paralog.candidates$Target, "_|_", sample, "_|_copy",
        sprintf("%02d", paralog.candidates$Candidate_rank))
      paralog.file = file.path(paralog.directory, paste0(sample, ".fa"))
      temp.paralog = tempfile(paste0(sample, "-"), tmpdir = paralog.directory)
      PhyloProcessR::writeFasta(
        sequences = as.list(as.character(paralog.seqs)),
        names = names(paralog.seqs), file = temp.paralog,
        nbchar = 1000000, as.string = TRUE)
      if (!file.rename(temp.paralog, paralog.file)) {
        unlink(temp.paralog)
        stop("Could not publish paralog FASTA for ", sample)
      }
    }

    print(paste0(sample, " target matching complete. ", nrow(selected),
                 " targets retained; ", length(unique(ambiguous.targets)),
                 " targets had competitive copies."))

    return(data.frame(
      Sample = sample,
      InputContigs = length(all.data),
      TargetsMatched = length(unique(qualifying$Target)),
      AnnotatedTargets = nrow(selected),
      AmbiguousTargets = length(unique(ambiguous.targets)),
      ParalogCandidates = nrow(paralog.candidates),
      MeanPident = if (nrow(qualifying) == 0) NA_real_ else
        round(mean(qualifying$Weighted_identity), 2),
      MeanBitscore = if (nrow(qualifying) == 0) NA_real_ else
        round(mean(qualifying$Total_bitscore), 1),
      MaxBitscore = if (nrow(qualifying) == 0) NA_real_ else
        round(max(qualifying$Total_bitscore), 1),
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    print(paste0(file.names[i], " failed: ", conditionMessage(e)))
    return("failed")
  })
  }, mc.cores = threads) # end i loop

  # A warning raised in a forked child never reaches the parent, so a failed
  # sample must be counted here.
  fail.count = sum(vapply(results, function(x) identical(x, "failed"), logical(1)))
  if (fail.count != 0) {
    print(paste0(fail.count, " of ", length(file.names),
                 " samples failed. See the messages above."))
  }

  ########################################################################
  # Write cross-sample summary log
  ########################################################################

  summary.df = do.call(rbind, results[vapply(results, is.data.frame, logical(1))])
  if (!is.null(summary.df) && nrow(summary.df) > 0) {
    out.csv = "logs/annotateTargets_summary.csv"
    if (file.exists(out.csv) && !overwrite) {
      existing = read.csv(out.csv, stringsAsFactors = FALSE)
      existing = existing[!existing$Sample %in% summary.df$Sample, ]
      summary.df = rbind(existing, summary.df)
    }
    write.csv(summary.df, file = out.csv, row.names = FALSE)
  }

  if (fail.count != 0) {
    stop(fail.count, " of ", length(file.names),
         " samples failed during annotateTargets; completed outputs were retained for resume.")
  }

  ########################################################################
  # Output a single file for alignment
  ########################################################################

  #gets lists of directories and files with sample names
  file.names = list.files(assembly.directory, pattern = "\\.fa$", full.names = FALSE)
  samples = gsub(".fa$", "", file.names)

  header.data = c("Sample", "startContigs", "annotatedContigs", "minLen", "maxLen", "meanLen")
  save.data = data.table::data.table(matrix(as.double(0), nrow = length(samples), ncol = length(header.data)))
  data.table::setnames(save.data, header.data)
  save.data[, Sample:=as.character(samples)]

  #Cycles through each assembly run and assesses each
  save.contigs = Biostrings::DNAStringSet()
  for (i in seq_along(samples)){

    #Gets length of contigs
    og.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", samples[i], ".fa"))
    data.table::set(save.data, i = i, j = match("startContigs", header.data),
                    value = length(og.contigs))
    out.fa = paste0(output.directory, "/", samples[i], ".fa")
    if (!file.exists(out.fa)) {
      warning(samples[i],
              ": no annotated output file found; no targets passed the filters.")
      next
    }
    cd.contigs = Biostrings::readDNAStringSet(out.fa)

    #Gets the saved matching targets
    data.table::set(save.data, i = i, j = match("Sample", header.data),
                    value = samples[i])
    data.table::set(save.data, i = i, j = match("annotatedContigs", header.data),
                    value = length(cd.contigs))

    contig.widths = Biostrings::width(cd.contigs)
    data.table::set(save.data, i = i, j = match("minLen", header.data),
                    value = min(contig.widths))
    data.table::set(save.data, i = i, j = match("maxLen", header.data),
                    value = max(contig.widths))
    data.table::set(save.data, i = i, j = match("meanLen", header.data),
                    value = mean(contig.widths))

    names(cd.contigs) = paste0(sub("_\\|_.*$", "", names(cd.contigs)),
                               "_|_", samples[i])
    save.contigs = append(save.contigs, cd.contigs)

  }#End loop for things

  #Finds probes that match to two or more contigs
  final.loci = as.list(as.character(save.contigs))
  combined.file = paste0(alignment.contig.name, "_to-align.fa")
  dir.create(dirname(combined.file), recursive = TRUE, showWarnings = FALSE)
  combined.temp = tempfile("annotated-targets-", tmpdir = dirname(combined.file))
  PhyloProcessR::writeFasta(sequences = final.loci, names = names(final.loci),
             combined.temp, nbchar = 1000000, as.string = T)
  if (!file.rename(combined.temp, combined.file)) {
    unlink(combined.temp)
    stop("Could not publish the combined annotation FASTA.")
  }

  #Saves combined, final dataset
  write.csv(save.data, file = "logs/annotation_sample_summary.csv", row.names = F)

} #End function


#### ** IDEA
#### Make table of contig matches to tagets, coordinates in each and length

#### END SCRIPT
