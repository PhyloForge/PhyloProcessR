#' @title pairwiseDistanceTarget
#'
#' @description Function for trimming out divergent sample sequence segments in an alignment
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param target a single sequence or reference aligned to alignment (or the name in the alignment) to pairwise compare each sample
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


#Finds the alignemnt pairwise distance from the target
pairwiseDistanceTarget = function(alignment = NULL,
                                  target = NULL) {

  #Alignment should be DNAStringSet
  #Debug
  # alignment = alignment
  #target = con.seq

  #Sets up reference
  if (class(target) == "DNAStringSet"){
    new.align = append(alignment, target)
    target.name = names(target)
  }
  if (class(target) == "character"){
    new.align  = alignment
    target.name = target
  }

  #Convert alignment to easier to work with matrix
  temp.align = strsplit(as.character(new.align), "")
  mat.align = lapply(temp.align, tolower)
  m.align = as.matrix(ape::as.DNAbin(mat.align))

  #Filters out weirdly divergent sequences
  n.align = as.character(m.align)
  n.align[n.align == "n"] = "-"
  n.align[n.align == "?"] = "-"
  n.align[is.na(n.align) == T] = "-"
  ref = n.align[rownames(n.align) == target.name,]

  #Processes each sample in the alignment
  summary.data = c()
  all.pars = c()
  all.over = c()
  for (z in 1:nrow(n.align)) {
    #Site counter
    pars = 0
    overlap = 0
    tar = n.align[z,]
    combined = matrix(NA_character_, ncol = max(length(ref), length(tar)), nrow =2)
    combined[1,] = ref
    combined[2,] = tar
    for (k in 1:ncol(combined)) {
      #Pulls out column of data
      seq.col = vector("character", length = nrow(combined))
      seq.col = combined[,k]
      #not equal to -
      f.char = seq.col[seq.col != '-']
      #don't count missing seq
      if (length(f.char) <= 1) { next }

      if (length(f.char) >= 2){
        overlap<-overlap+1
        if (f.char[1] != f.char [2]) { pars = pars + 1 }
      }#end if
    }#ends informative sites loop
    all.pars = append(all.pars, pars)
    all.over = append(all.over, overlap)
  }# ends seq loop

  #Summarizes and returns data
  summary.data = all.pars/all.over
  summary.data[is.nan(summary.data)] = 0
  names(summary.data) = names(new.align)

  #Returns the summary data
  return(summary.data)

}#end function
