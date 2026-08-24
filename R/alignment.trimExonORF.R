#' @title trimExonORF
#'
#' @description Optimizes and cleans exon/CDS multiple sequence alignments for protein-coding analysis. Evaluates all 6 reading frames (+1, +2, +3, -1, -2, -3) to detect the consensus open reading frame, optionally flips reverse-strand exons to the coding (+) strand, trims leading/trailing partial codons to preserve exact triplet boundaries (modulo 3 equals zero), and prunes or masks samples containing premature internal stop codons. Non-coding or UCE alignments can be automatically excluded.
#'
#' @param alignment a DNAStringSet containing the aligned exon sequences to process
#'
#' @param locus.name character name or identifier of the locus (used to detect and exclude UCE markers)
#'
#' @param auto.shift.frame if TRUE, trim leading offset (1-2 bp) and trailing partial codons so column 1 is codon position 1 and total length is divisible by 3 (default: TRUE)
#'
#' @param auto.flip.reverse if TRUE, automatically reverse complement the alignment if the optimal open reading frame is on the minus strand (default: TRUE)
#'
#' @param stop.codon.action action to take on samples with premature internal stop codons: "remove.sample" (default), "mask.codon" (replace stop with "---"), or "keep"
#'
#' @param genetic.code genetic code table to use for translation: "standard" (default), "vertebrate.mitochondrial", or "invertebrate.mitochondrial"
#'
#' @param exclude.uce if TRUE, alignments whose locus.name begins with "uce" or contains "noncoding" are returned unmodified (default: TRUE)
#'
#' @return a DNAStringSet with reading frame optimized and stop-codon-containing samples pruned or masked
#'
#' @export

trimExonORF = function(alignment = NULL,
                       locus.name = NULL,
                       auto.shift.frame = TRUE,
                       auto.flip.reverse = TRUE,
                       stop.codon.action = c("remove.sample", "mask.codon", "keep"),
                       genetic.code = c("standard", "vertebrate.mitochondrial", "invertebrate.mitochondrial"),
                       exclude.uce = TRUE) {

  #Debug section
  # alignment = align
  # locus.name = "exon_001"
  # auto.shift.frame = TRUE
  # auto.flip.reverse = TRUE
  # stop.codon.action = "remove.sample"
  # genetic.code = "standard"
  # exclude.uce = TRUE

  if (length(alignment) <= 2){ return(alignment) }

  stop.codon.action = match.arg(stop.codon.action)
  genetic.code = match.arg(genetic.code)

  #Check for UCE exclusion
  if (exclude.uce == TRUE && !is.null(locus.name)){
    lower_name = tolower(locus.name)
    if (grepl("^uce", lower_name) || grepl("noncoding", lower_name) || grepl("intron", lower_name)){
      return(alignment)
    }
  }

  #Helper function for reverse complement
  rev_comp_dna = function(seq_str){
    chars = strsplit(seq_str, "")[[1]]
    comp_map = c("A"="T", "T"="A", "U"="A", "C"="G", "G"="C",
                 "a"="t", "t"="a", "u"="a", "c"="g", "g"="c",
                 "-"="-", "?"="?", "N"="N", "n"="n",
                 "R"="Y", "Y"="R", "S"="S", "W"="W", "K"="M", "M"="K")
    rev_chars = rev(chars)
    mapped = comp_map[rev_chars]
    mapped[is.na(mapped)] = "N"
    return(paste0(mapped, collapse = ""))
  }

  #Helper for single codon translation
  translate_triplet = function(trip_str, code_type){
    if (nchar(trip_str) < 3){ return("-") }
    up = toupper(trip_str)
    if (grepl("-", up)){ return("-") }
    if (grepl("[?N]", up)){ return("X") }

    if (code_type == "standard"){
      if (up %in% c("TAA", "TAG", "TGA")){ return("*") }
    } else if (code_type == "vertebrate.mitochondrial"){
      if (up %in% c("TAA", "TAG", "AGA", "AGG")){ return("*") }
      if (up == "TGA"){ return("W") }
      if (up == "ATA"){ return("M") }
    } else if (code_type == "invertebrate.mitochondrial"){
      if (up %in% c("TAA", "TAG")){ return("*") }
      if (up == "TGA"){ return("W") }
      if (up == "ATA"){ return("M") }
      if (up %in% c("AGA", "AGG")){ return("S") }
    }
    #Standard amino acid mappings
    std_table = c(
      "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L", "CTT"="L", "CTC"="L", "CTA"="L", "CTG"="L",
      "ATT"="I", "ATC"="I", "ATA"="I", "ATG"="M", "GTT"="V", "GTC"="V", "GTA"="V", "GTG"="V",
      "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S", "AGT"="S", "AGC"="S", "CCT"="P", "CCC"="P",
      "CCA"="P", "CCG"="P", "ACT"="T", "ACC"="T", "ACA"="T", "ACG"="T", "GCT"="A", "GCC"="A",
      "GCA"="A", "GCG"="A", "TAT"="Y", "TAC"="Y", "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
      "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K", "GAT"="D", "GAC"="D", "GAA"="E", "GAG"="E",
      "TGT"="C", "TGC"="C", "TGG"="W", "CGT"="R", "CGC"="R", "CGA"="R", "CGG"="R", "AGA"="R",
      "AGG"="R", "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G"
    )
    aa = std_table[up]
    if (is.na(aa)){ return("X") }
    return(aa)
  }

  raw_seqs = as.character(alignment)
  rev_seqs = unname(sapply(raw_seqs, rev_comp_dna))

  #Find optimal reading frame among 6 possible frames
  best_frame = 1
  best_is_rev = FALSE
  best_offset = 0
  best_score = -Inf

  for (is_rev in c(FALSE, TRUE)){
    test_seqs = if (is_rev) rev_seqs else raw_seqs

    for (offset in 0:2){
      frame_num = if (is_rev) -(offset + 1) else (offset + 1)
      total_clean_codons = 0
      start_codons = 0
      terminal_stops = 0
      stop_taxa = 0

      for (seq in test_seqs){
        seq_len = nchar(seq)
        if (seq_len < offset + 3){ next }
        codon_count = floor((seq_len - offset) / 3)
        has_stop = FALSE

        for (c_idx in 0:(codon_count - 1)){
          start_pos = offset + c_idx * 3 + 1
          trip = substr(seq, start_pos, start_pos + 2)
          aa = translate_triplet(trip, genetic.code)

          if (c_idx == 0 && toupper(trip) == "ATG"){ start_codons = start_codons + 1 }
          if (c_idx == (codon_count - 1) && aa == "*"){ terminal_stops = terminal_stops + 1 }
          if (aa == "*" && c_idx < (codon_count - 1)){ has_stop = TRUE }
        }

        if (has_stop == TRUE){
          stop_taxa = stop_taxa + 1
        } else {
          total_clean_codons = total_clean_codons + codon_count
        }
      }#end seq loop

      score = (total_clean_codons * 100) + (start_codons * 250) + (terminal_stops * 150) - (stop_taxa * 80)
      if (score > best_score){
        best_score = score
        best_frame = frame_num
        best_is_rev = is_rev
        best_offset = offset
      }
    }#end offset loop
  }#end is_rev loop

  #Apply strand flipping if needed
  curr_seqs = if (best_is_rev == TRUE && auto.flip.reverse == TRUE) rev_seqs else raw_seqs
  curr_names = names(alignment)

  #Apply frame shift and codon triplet length trimming
  if (auto.shift.frame == TRUE){
    old_len = max(nchar(curr_seqs))
    trailing_rem = if (old_len >= best_offset) (old_len - best_offset) %% 3 else 0
    new_len = max(0, old_len - best_offset - trailing_rem)
    start_pos = best_offset + 1
    end_pos = start_pos + new_len - 1

    if (new_len > 0){
      curr_seqs = substr(curr_seqs, start_pos, end_pos)
    } else {
      return(alignment)
    }
  }

  #Sample-level stop codon filtering or masking
  if (stop.codon.action != "keep"){
    align_len = nchar(curr_seqs[1])
    codon_count = floor(align_len / 3)

    kept_indices = c()
    final_seqs = c()

    for (s_idx in seq_along(curr_seqs)){
      seq_str = curr_seqs[s_idx]
      seq_chars = strsplit(seq_str, "")[[1]]
      has_internal_stop = FALSE

      if (codon_count > 0){
        for (c_idx in 0:(codon_count - 1)){
          start_pos = c_idx * 3 + 1
          trip = substr(seq_str, start_pos, start_pos + 2)
          aa = translate_triplet(trip, genetic.code)

          if (aa == "*" && c_idx < (codon_count - 1)){
            has_internal_stop = TRUE
            if (stop.codon.action == "mask.codon"){
              seq_chars[start_pos:(start_pos + 2)] = "-"
            }
          }
        }
      }

      if (has_internal_stop == TRUE && stop.codon.action == "remove.sample"){
        #Drop sample
      } else {
        kept_indices = c(kept_indices, s_idx)
        final_seqs = c(final_seqs, paste0(seq_chars, collapse = ""))
      }
    }#end s_idx loop

    curr_seqs = final_seqs
    curr_names = curr_names[kept_indices]
  }

  if (length(curr_seqs) == 0){ return(Biostrings::DNAStringSet()) }

  output_align = Biostrings::DNAStringSet(curr_seqs)
  names(output_align) = curr_names

  return(output_align)

}#end trimExonORF
