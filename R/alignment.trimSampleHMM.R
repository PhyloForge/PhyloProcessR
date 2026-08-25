#' @title trimSampleHMM
#'
#' @description Masks aberrant, misaligned, or chimeric sequence segments within individual samples in an alignment using Leave-One-Out Profile HMM posterior match decoding (similar to TAPIR / HMMCleaner). For each column in the alignment, a position-specific emission profile is computed from all other samples (with Dirichlet smoothing), and the posterior match confidence is evaluated for every residue. A rolling average of the posterior is computed over a window to prevent random sequence matches from disrupting the detection of contiguous bad segments. Windows whose average confidence falls below `min.posterior` are replaced with gaps. Alignments with two or fewer sequences are returned unmodified.
#'
#' @param alignment a DNAStringSet containing the aligned sequences to clean
#'
#' @param min.posterior numeric threshold (0-1): residues with a smoothed Leave-One-Out Profile HMM posterior match confidence below this value are considered aberrant (default: 0.45)
#'
#' @param min.segment.length integer: the window size for smoothing and the minimum length of aberrant segments targeted for masking (default: 8)
#'
#' @param min.island.length integer: minimum span of unmasked bases required to be kept between masked segments (default: 20)
#'
#' @param mask.char character to replace aberrant bases with; defaults to "-"
#'
#' @return a DNAStringSet with aberrant sequence segments masked as gaps
#'
#' @export

trimSampleHMM = function(alignment = NULL,
                         min.posterior = 0.45,
                         min.segment.length = 8,
                         min.island.length = 20,
                         mask.char = "-") {

  if (length(alignment) <= 2){ return(alignment) }

  #Converts DNAStringSet to uppercase character matrix
  seq_chars = strsplit(as.character(alignment), "")
  char_matrix = do.call(rbind, seq_chars)
  upper_matrix = toupper(char_matrix)

  num_taxa = nrow(char_matrix)
  num_cols = ncol(char_matrix)

  if (num_cols == 0){ return(alignment) }

  alpha = 0.5
  null_prob = 0.25

  #Precomputes column-wise observed base counts (A, C, G, T)
  counts_A = colSums(upper_matrix == "A")
  counts_C = colSums(upper_matrix == "C")
  counts_G = colSums(upper_matrix == "G")
  counts_T = colSums(upper_matrix == "T" | upper_matrix == "U")
  col_totals = counts_A + counts_C + counts_G + counts_T

  out_matrix = char_matrix
  half_w = floor(min.segment.length / 2)

  #Scores each sample using Leave-One-Out Jackknife Profile HMM
  for (i in 1:num_taxa) {
    seq_i = upper_matrix[i, ]
    confidences = rep(1.0, num_cols)

    is_A = seq_i == "A"
    is_C = seq_i == "C"
    is_G = seq_i == "G"
    is_T = seq_i == "T" | seq_i == "U"
    is_valid = is_A | is_C | is_G | is_T

    #Leave-one-out counts: exclude sample i from its own profile
    count_without_i = numeric(num_cols)
    count_without_i[is_A] = pmax(0, counts_A[is_A] - 1)
    count_without_i[is_C] = pmax(0, counts_C[is_C] - 1)
    count_without_i[is_G] = pmax(0, counts_G[is_G] - 1)
    count_without_i[is_T] = pmax(0, counts_T[is_T] - 1)

    total_without_i = pmax(0, col_totals - 1)
    denom = total_without_i + alpha
    p_emit = (count_without_i + alpha * null_prob) / denom

    #Posterior match confidence
    confidences[is_valid] = p_emit[is_valid] / (p_emit[is_valid] + null_prob)

    # Calculate rolling mean over window (O(N) using cumsum)
    cum_conf = c(0, cumsum(confidences))
    start_cols = pmax(1, 1:num_cols - half_w)
    end_cols = pmin(num_cols, 1:num_cols + half_w)
    
    smoothed = (cum_conf[end_cols + 1] - cum_conf[start_cols]) / (end_cols - start_cols + 1)
    
    # Identify windows that drop below threshold
    bad_centers = which(smoothed < min.posterior)
    if (length(bad_centers) > 0) {
      out_matrix[i, bad_centers] = mask.char
    }
    
    if (min.island.length > 0) {
      is_masked = out_matrix[i, ] == mask.char
      runs = rle(is_masked)
      
      if (length(runs$lengths) > 1) {
        run_starts = c(1, cumsum(runs$lengths)[-length(runs$lengths)] + 1)
        for (r in 1:length(runs$lengths)) {
          # if it's an island of FALSE (not masked) and shorter than island length
          if (!runs$values[r] && runs$lengths[r] < min.island.length) {
            start_idx = run_starts[r]
            end_idx = run_starts[r] + runs$lengths[r] - 1
            out_matrix[i, start_idx:end_idx] = mask.char
          }
        }
      }
    }
  }#end i loop

  #Converts matrix back to DNAStringSet preserving names
  res_strings = apply(out_matrix, 1, paste0, collapse = "")
  output_align = Biostrings::DNAStringSet(res_strings)
  names(output_align) = names(alignment)

  return(output_align)

}#end trimSampleHMM
