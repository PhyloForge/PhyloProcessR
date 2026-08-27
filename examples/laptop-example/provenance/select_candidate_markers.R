args <- commandArgs(trailingOnly = TRUE)
counts <- read.delim(args[[1]], stringsAsFactors = FALSE)
samples <- sort(unique(counts$sample))

wide <- reshape(
  counts,
  idvar = "marker",
  timevar = "sample",
  direction = "wide"
)
count_cols <- paste0("mapped_alignments.", samples)
for (column in count_cols) {
  if (!column %in% names(wide)) wide[[column]] <- 0L
  wide[[column]][is.na(wide[[column]])] <- 0L
}

wide$samples_supported <- rowSums(wide[count_cols] > 0)
wide$minimum_count <- apply(wide[count_cols], 1, min)
wide$total_count <- rowSums(wide[count_cols])
wide$gene <- sub("-ex.*$", "", sub("^[^_]+_", "", wide$marker))
wide <- wide[wide$samples_supported == length(samples), ]
wide <- wide[order(-wide$minimum_count, -wide$total_count, wide$marker), ]
wide <- wide[!duplicated(wide$gene), ]

forced <- c("ranoidea-01644_kif2a-ex2", "ranoidea-05399_ckap5-ex3")
selected <- unique(c(forced, wide$marker))
selected <- selected[seq_len(min(40, length(selected)))]

write.table(
  wide,
  args[[2]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
writeLines(selected, args[[3]])

cat("samples:", length(samples), "\n")
cat("markers present in all samples:", nrow(wide), "gene-unique\n")
cat("selected candidate markers:", length(selected), "\n")
cat("weakest selected minimum count:", min(wide$minimum_count[match(setdiff(selected, forced), wide$marker)], na.rm = TRUE), "\n")
