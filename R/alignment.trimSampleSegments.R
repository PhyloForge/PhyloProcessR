#' @title trimSampleSegments
#'
#' @description Function for trimming out divergent sample sequence segments in an alignment
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param slice.size.bp size of each sliding window to calculate divergence
#'
#' @param threshold pair-wise distance threshold from consensus to trim out segment
#'
#' @return returns DNAStringSet of trimmed alignment
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export

#Slice function
trimSampleSegments = function(alignment = NULL,
                              reference.name = NULL,
                              slice.size.bp = 100,
                              threshold = 0.45){

  #makes consensus sequence for comparison
  #input.align<-align
  #reference.name = names(add.ref)
  #slice.size.bp = 50
  #threshold = 0.4

  if (length(alignment) <= 3){ return(alignment) }

  if (is.null(reference.name) != T){
    comb.align = alignment
  } else {
    input.con = PhyloCap::makeConsensus(alignment = alignment,
                              method = "majority",
                              remove.gaps = F)

    names(input.con) = "Reference_Locus"
    comb.align = append(alignment, input.con)
    reference.name = "Reference_Locus"
  }#end else

  #Gets slice information ready
  slice.no = ceiling(max(Biostrings::width(alignment))/slice.size.bp)
  slice.start = 1
  slice.end = slice.size.bp

  #checks to see if its out of bounds
  if (slice.end > max(Biostrings::width(alignment))){
    slice.end = max(Biostrings::width(alignment))
  }#end if check

  #Loops through each slice
  output.align = Biostrings::DNAStringSet()
  for (x in 1:slice.no){

    #Slice alignment into number of slices
    sliced.align = subseq(comb.align, start = slice.start, end = slice.end)

    #Checks for badly aligned sequences
    bad.align = PhyloCap::pairwiseDistanceTarget(sliced.align, reference.name)

    #Remove bad sequence chunks
    rem.seqs = bad.align[bad.align >= threshold]
    good.align = sliced.align[!names(sliced.align) %in% names(rem.seqs)]

    #Makes replacement gap seqs for the bad ones
    blank.align = DNAStringSet()
    if (length(rem.seqs) != 0){
      for (y in 1:length(rem.seqs)){
        blank.align = append(blank.align, Biostrings::DNAStringSet(paste0(rep("-", slice.end-slice.start+1), collapse = "")) )
      }
      names(blank.align) = names(rem.seqs)
    }#end rem seqs if

    #Saves the slices and cats
    save.slice = append(good.align, blank.align)
    save.slice = save.slice[order(names(save.slice))]
    save.names = names(save.slice)
    output.align = Biostrings::DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
    names(output.align) = save.names

    #Gets new start and stop
    slice.start = slice.start+slice.size.bp
    slice.end = slice.end+slice.size.bp

    #checks to see if its out of bounds
    if (slice.end > max(Biostrings::width(alignment))){
      slice.end = max(Biostrings::width(alignment))
    }#end oversized if

    if (slice.end <= slice.start){ break }

    if (slice.end-slice.start < threshold){
      save.slice = Biostrings::subseq(comb.align, start = slice.start, end = slice.end)
      save.slice = save.slice[order(names(save.slice))]
      save.names = names(save.slice)
      output.align = Biostrings::DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
      names(output.align) = save.names
      break
    }#end if

  }#end x loop

  #Removes reference
  if (reference.name == "Reference_Locus"){
    output.align = output.align[names(output.align) != "Reference_Locus"]
  }#end removal

  #removes gap only taxa
  str.splitted = strsplit(as.character(output.align), "")
  x.align = as.matrix(ape::as.DNAbin(str.splitted) )
  len.temp = as.character(as.list(x.align))
  len.loci = lapply(len.temp, function (x) x[x != "-"])
  spp.len = unlist(lapply(len.loci, function (x) length(x)))
  spp.rem = spp.len[spp.len <= slice.size.bp]
  return.align = output.align[!names(output.align) %in% names(spp.rem)]

  return(return.align)

}#end FUNCTION

