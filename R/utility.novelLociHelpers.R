# Find bases covered in enough samples before joining nearby shared intervals.
# Column 4 records the minimum support across the retained covered intervals;
# gaps introduced by max.merge.distance do not have this coverage guarantee.
.sharedCoveredRegions = function(bed.files, shared.bed, min.samples,
                                  min.region.length, max.merge.distance,
                                  bedtools.command, quiet = TRUE) {

  non.empty = bed.files[file.size(bed.files) > 0]
  if (length(non.empty) < min.samples || length(non.empty) == 0) {
    writeLines(character(), shared.bed)
    return(invisible(NULL))
  }

  # multiinter requires the same chromosome sort order in every sample.
  sorted.files = paste0(non.empty, ".sorted")
  on.exit(unlink(sorted.files), add = TRUE)
  for (i in seq_along(non.empty)) {
    .runCommand(paste0("LC_ALL=C sort -k1,1 -k2,2n ", shQuote(non.empty[i]),
                        " > ", shQuote(sorted.files[i])),
                quiet = quiet, keep.stdout = TRUE, task = "coverage BED sorting")
  }

  if (length(sorted.files) == 1) {
    # multiinter requires at least two input files.
    intersect.command = paste0("awk 'BEGIN{OFS=\"\\t\"}{print $1,$2,$3,1}' ",
                                 shQuote(sorted.files))
  } else {
    intersect.command = paste0(bedtools.command, " multiinter -i ",
                                 paste(shQuote(sorted.files), collapse = " "))
  }
  temp.bed = paste0(shared.bed, ".tmp")
  on.exit(unlink(temp.bed), add = TRUE)
  .runPipeline(paste0(intersect.command,
                       " | awk 'BEGIN{OFS=\"\\t\"} $4 >= ", min.samples,
                       " {print $1,$2,$3,$4}' | ", bedtools.command,
                       " merge -i stdin -d ", max.merge.distance, " -c 4 -o min",
                       " | awk '($3-$2) >= ", min.region.length, "' > ", shQuote(temp.bed)),
               quiet = quiet, keep.stdout = TRUE, task = "shared region selection")
  if (!file.rename(temp.bed, shared.bed)) stop("Cannot save shared region BED.")
  invisible(NULL)
}#end .sharedCoveredRegions
