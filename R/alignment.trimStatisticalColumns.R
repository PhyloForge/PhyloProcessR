#' @title trimStatisticalColumns
#'
#' @description Trims alignment columns based on statistical properties such as pairwise similarity (trimAl style), conserved block structures (Gblocks style), or Shannon information entropy.
#'
#' @param alignment a DNAStringSet containing the aligned sequences.
#' @param method character: "trimal", "gblocks", or "entropy".
#' @param similarity.threshold numeric: threshold for column similarity (used in trimal custom heuristic or gblocks).
#' @param window.size integer: sliding window size for smoothing similarity scores (trimAl).
#' @param heuristic character: "custom", "gappyout", "strict", or "strictplus" (trimAl heuristics).
#' @param min.block.length integer: minimum number of contiguous conserved columns to form a block (Gblocks).
#' @param max.nonconserved integer: maximum number of contiguous non-conserved columns allowed within a block (Gblocks).
#' @param gap.treatment character: "none", "half", or "all" determining how gaps affect conserved status (Gblocks).
#' @param max.entropy numeric: maximum allowed Shannon entropy in bits (Entropy method).
#'
#' @return a DNAStringSet with the filtered columns removed.
#'
#' @export

trimStatisticalColumns = function(alignment = NULL,
                                  method = c("trimal", "gblocks", "entropy"),
                                  similarity.threshold = 0.35,
                                  window.size = 3,
                                  heuristic = c("custom", "gappyout", "strict", "strictplus"),
                                  min.block.length = 5,
                                  max.nonconserved = 4,
                                  gap.treatment = c("half", "none", "all"),
                                  max.entropy = 1.5) {
                                  
  method = match.arg(method)
  heuristic = match.arg(heuristic)
  gap.treatment = match.arg(gap.treatment)
  
  if (length(alignment) == 0) return(alignment)
  
  seq_chars = strsplit(as.character(alignment), "")
  char_matrix = do.call(rbind, seq_chars)
  char_matrix = toupper(char_matrix)
  
  num_taxa = nrow(char_matrix)
  num_cols = ncol(char_matrix)
  
  if (num_cols == 0) return(alignment)
  
  # Helper to compute column gap fraction
  compute_gap_fraction = function(mat) {
    colSums(mat == "-" | mat == "?" | mat == "N") / nrow(mat)
  }
  
  gap_fracs = compute_gap_fraction(char_matrix)
  
  # Helper to compute pairwise similarity
  compute_similarity = function(mat) {
    apply(mat, 2, function(col) {
      valid = col[col != "-" & col != "?" & col != "N"]
      n = length(valid)
      if (n <= 1) return(0)
      counts = table(valid)
      matches = sum(counts * (counts - 1) / 2)
      total = n * (n - 1) / 2
      matches / total
    })
  }
  
  # Helper to apply sliding window smoothing
  apply_window = function(scores, w) {
    if (w <= 1) return(scores)
    half = floor(w / 2)
    sapply(1:length(scores), function(i) {
      start = max(1, i - half)
      end = min(length(scores), i + half)
      mean(scores[start:end])
    })
  }
  
  keep_cols = rep(TRUE, num_cols)
  
  if (method == "trimal") {
    raw_scores = compute_similarity(char_matrix)
    smoothed = apply_window(raw_scores, window.size)
    
    avg_score = mean(smoothed)
    
    eff_thresh = switch(heuristic,
                        "custom" = similarity.threshold,
                        "gappyout" = max(0.15, min(0.60, avg_score * 0.75)),
                        "strict" = max(0.25, min(0.70, avg_score * 0.90)),
                        "strictplus" = max(0.35, min(0.85, avg_score * 1.10)))
                        
    for (i in 1:num_cols) {
      is_dropped = FALSE
      if (heuristic == "gappyout") {
        is_dropped = (gap_fracs[i] > 0.60) || (smoothed[i] < eff_thresh)
      } else if (heuristic %in% c("strict", "strictplus")) {
        is_dropped = (gap_fracs[i] > 0.40) || (smoothed[i] < eff_thresh)
      } else {
        is_dropped = (smoothed[i] < eff_thresh)
      }
      keep_cols[i] = !is_dropped
    }
    
  } else if (method == "gblocks") {
    raw_scores = compute_similarity(char_matrix)
    
    is_conserved = rep(FALSE, num_cols)
    for (i in 1:num_cols) {
      gap_pass = switch(gap.treatment,
                        "none" = (gap_fracs[i] == 0),
                        "half" = (gap_fracs[i] <= 0.5),
                        "all" = TRUE)
      if (gap_pass && raw_scores[i] >= similarity.threshold) {
        is_conserved[i] = TRUE
      }
    }
    
    # Extract blocks
    in_block = rep(FALSE, num_cols)
    current_block = c()
    non_conserved_count = 0
    
    for (i in 1:num_cols) {
      if (is_conserved[i]) {
        current_block = c(current_block, i)
        non_conserved_count = 0
      } else {
        if (length(current_block) > 0) {
          non_conserved_count = non_conserved_count + 1
          current_block = c(current_block, i)
          if (non_conserved_count > max.nonconserved) {
            # Block broken. Check if valid before discarding non-conserved tail
            valid_len = length(current_block) - non_conserved_count
            if (valid_len >= min.block.length) {
              in_block[current_block[1:(valid_len)]] = TRUE
            }
            current_block = c()
            non_conserved_count = 0
          }
        }
      }
    }
    
    if (length(current_block) > 0) {
      valid_len = length(current_block) - non_conserved_count
      if (valid_len >= min.block.length) {
        in_block[current_block[1:(valid_len)]] = TRUE
      }
    }
    
    keep_cols = in_block
    
  } else if (method == "entropy") {
    compute_entropy = function(mat) {
      apply(mat, 2, function(col) {
        valid = col[col %in% c("A", "C", "G", "T", "U")]
        n = length(valid)
        if (n <= 1) return(0)
        counts = table(valid)
        p = counts / n
        -sum(p * log2(p))
      })
    }
    
    entropies = compute_entropy(char_matrix)
    keep_cols = entropies <= max.entropy
  }
  
  if (!any(keep_cols)) {
    # If all columns removed, return empty DNAStringSet but preserve taxa names
    empty_strings = rep("", num_taxa)
    output = Biostrings::DNAStringSet(empty_strings)
    names(output) = names(alignment)
    return(output)
  }
  
  filtered_matrix = char_matrix[, keep_cols, drop = FALSE]
  res_strings = apply(filtered_matrix, 1, paste0, collapse = "")
  output = Biostrings::DNAStringSet(res_strings)
  names(output) = names(alignment)
  
  return(output)
}
