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
                              slice.size.bp = 100,
                              threshold = 0.45){

  #makes consensus sequence for comparison
  #input.align<-align

  if (length(alignment) <= 2){ return(alignment) }

  input.con = makeConsensus(alignment = alignment,
                            method = "majority")

  names(input.con) = "Reference_Locus"

  # #B Slice up alignment [combine with prev]
  # if (max(width(trimal.align)) >= 100){
  #   red.align<-slice.trim(trimal.align, slice.size.bp = 80, threshold = 0.40)
  #   if (length(red.align) == 0){ next }
  # }else {red.align = trimal.align }

  comb.align = append(input.align, input.con)

  #Gets slice information ready
  slice.no = ceiling(max(width(input.align))/slice.size.bp)
  slice.start = 1
  slice.end = slice.size.bp

  #checks to see if its out of bounds
  if (slice.end > max(width(input.align))){
    slice.end<-max(width(input.align))
  }#end if check
  output.align = DNAStringSet()
  for (x in 1:slice.no){

    #Slice alignment into number of slices
    sliced.align = subseq(comb.align, start = slice.start, end = slice.end)
    #Checks for badly aligned sequences
    bad.align = pairwise.inf.sites(sliced.align, "Reference_Locus")
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
    slice.start = slice.start+100
    slice.end = slice.end+100
    #checks to see if its out of bounds
    if (slice.end > max(width(input.align))){
      slice.end<-max(width(input.align))
      if (slice.end-slice.start <= 25){ break } else {
        save.slice<-subseq(comb.align, start = slice.start, end = slice.end)
        save.slice<-save.slice[order(names(save.slice))]
        save.names<-names(save.slice)
        output.align<-Biostrings::DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
        names(output.align)<-save.names
        break
      }
    }#end if
  }#end x loop

  #Removes reference
  output.align<-output.align[names(output.align) != "Reference_Locus"]
  #removes gap only taxa
  str.splitted<-strsplit(as.character(output.align), "")
  x.align<-as.matrix(as.DNAbin(str.splitted) )
  len.temp<-as.character(as.list(x.align))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= 20]
  return.align<-output.align[!names(output.align) %in% names(spp.rem)]

  return(return.align)

}#end FUNCTION

