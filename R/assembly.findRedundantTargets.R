#' @title findRedundantTargets
#'
#' @description Finds target markers that capture the same region of the genome.
#' A probe set made from several sources (for example UCE, AHE and BUSCO probes)
#' can have two targets for one region. In each sample, the function finds the
#' best contig of each target and the region of the contig that the target
#' matches. Two targets share a region when they have the same best contig and
#' their regions overlap by more than half of the shorter region. Two targets
#' are redundant when they share a region in enough samples. Redundant targets are put in groups, and each
#' group gets one representative target. \code{annotateTargets} reads the result
#' and gives the region to the representative in every sample, so the data of
#' one region stays in one alignment.
#'
#' @details The function reads the per-sample match logs that
#' \code{curateTargetContigs} writes (\code{<sample>_curation-matches.csv}). Run
#' it after the last curation step of workflow 2. Two targets that are genes of
#' one family usually have different best contigs, and two targets for
#' neighbouring exons match different regions of one contig, so neither is
#' grouped.
#'
#' @param log.directory path to the directory that contains the
#' \code{<sample>_curation-matches.csv} files. Default "logs/sample_logs".
#'
#' @param output.directory path to the directory for the two output tables.
#' Default "redundant-targets".
#'
#' @param min.samples minimum number of samples in which two targets must share
#' a region. Default 3.
#'
#' @param min.shared.proportion minimum proportion of the samples that have
#' both targets in which the two targets share a region. Default 0.5.
#'
#' @return Writes two files to \code{output.directory}:
#' \itemize{
#'   \item \code{redundant-targets.csv} -- one row per grouped target: the
#'     target, its representative, the group size, and the number of samples
#'     with the target. \code{annotateTargets} reads this file.
#'   \item \code{shared-contig-pairs.csv} -- one row per target pair that
#'     shares a region in at least one sample, with the counts and whether the
#'     pair was grouped. Use it to review pairs below the thresholds.
#' }
#' Returns the \code{redundant-targets.csv} table invisibly.
#'
#' @export

findRedundantTargets = function(log.directory = "logs/sample_logs",
                                output.directory = "redundant-targets",
                                min.samples = 3,
                                min.shared.proportion = 0.5) {

  # Parameter checks
  if (!dir.exists(log.directory)) stop("log.directory not found: ", log.directory)
  if (!is.numeric(min.samples) || length(min.samples) != 1 || min.samples < 1) {
    stop("min.samples must be one number of 1 or more.")
  }
  if (!is.numeric(min.shared.proportion) || length(min.shared.proportion) != 1 ||
      min.shared.proportion <= 0 || min.shared.proportion > 1) {
    stop("min.shared.proportion must be one number above 0 and up to 1.")
  }

  log.files = list.files(log.directory, pattern = "_curation-matches\\.csv$",
                         full.names = TRUE)
  if (length(log.files) == 0) {
    stop("No _curation-matches.csv files in ", log.directory,
         ". Run curateTargetContigs first.")
  }
  dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)

  # Step 1: find the best contig of each target in each sample, and the region
  # of that contig that the target matches. In the log, qName is the target and
  # tName is the contig.
  best.list = list()
  for (log.file in log.files) {
    hits = data.table::fread(log.file)
    if (nrow(hits) == 0) next
    if (!all(c("tStart", "tEnd") %in% names(hits))) {
      stop(log.file, " has no contig coordinates. Run curateTargetContigs again.")
    }
    scores = hits[, .(score = sum(bitscore),
                      start = min(pmin(tStart, tEnd)),
                      end = max(pmax(tStart, tEnd))),
                  by = .(Sample, qName, tName)]
    data.table::setorder(scores, qName, -score, tName)
    best.list[[log.file]] = scores[!duplicated(qName)]
  }
  best = data.table::rbindlist(best.list)
  data.table::setnames(best, c("qName", "tName"), c("Target", "Contig"))

  # Step 2: count, for each target pair, the samples in which the two targets
  # share a region: the same best contig, and an overlap of more than half of
  # the shorter region.
  pairs = merge(best, best, by = c("Sample", "Contig"),
                suffixes = c("_a", "_b"), allow.cartesian = TRUE)
  pairs = pairs[Target_a < Target_b]
  overlap = pmin(pairs$end_a, pairs$end_b) - pmax(pairs$start_a, pairs$start_b) + 1
  shorter = pmin(pairs$end_a - pairs$start_a + 1, pairs$end_b - pairs$start_b + 1)
  pairs = pairs[overlap > shorter / 2]
  shared = pairs[, .(Shared_samples = .N), by = .(Target_a, Target_b)]

  # Step 3: count the samples that have both targets, then mark the pairs that
  # pass both thresholds as redundant.
  samples.of = split(best$Sample, best$Target)
  shared$Samples_with_both = mapply(function(a, b) {
    length(intersect(samples.of[[a]], samples.of[[b]]))
  }, shared$Target_a, shared$Target_b)
  shared$Shared_proportion = shared$Shared_samples / shared$Samples_with_both
  shared$Grouped = shared$Shared_samples >= min.samples &
    shared$Shared_proportion >= min.shared.proportion
  data.table::setorder(shared, -Shared_samples, Target_a, Target_b)

  # Step 4: join the redundant pairs into groups. A chain of pairs (a-b, b-c)
  # gives one group (a, b, c). Each target starts in its own group, and the
  # loop moves both targets of a pair to the lower group name until nothing
  # changes.
  redundant = shared[Grouped == TRUE]
  grouped.targets = sort(unique(c(redundant$Target_a, redundant$Target_b)))
  group = stats::setNames(grouped.targets, grouped.targets)
  changed = TRUE
  while (changed) {
    changed = FALSE
    for (row in seq_len(nrow(redundant))) {
      a = redundant$Target_a[row]
      b = redundant$Target_b[row]
      lowest = min(group[[a]], group[[b]])
      if (group[[a]] != lowest || group[[b]] != lowest) {
        group[[a]] = lowest
        group[[b]] = lowest
        changed = TRUE
      }
    }
  }

  # Step 5: the representative of a group is the target found in the most
  # samples. When the counts are equal, the name that sorts first is used.
  result = data.table::data.table(
    Target = grouped.targets,
    Group = unname(group[grouped.targets]),
    Samples_with_target = lengths(samples.of[grouped.targets]))
  data.table::setorder(result, Group, -Samples_with_target, Target)
  result[, Representative := Target[1], by = Group]
  result[, Group_size := .N, by = Group]
  result = result[, .(Target, Representative, Group_size, Samples_with_target)]

  write.csv(result, file.path(output.directory, "redundant-targets.csv"),
            row.names = FALSE)
  write.csv(shared, file.path(output.directory, "shared-contig-pairs.csv"),
            row.names = FALSE)

  print(paste0(length(unique(result$Representative)), " groups of redundant targets (",
               nrow(result), " targets) from ", length(log.files), " samples."))
  return(invisible(result))
}
