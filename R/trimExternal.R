#' @title trimExternal
#'
#' @description Function for externally trimming at the edges in an alignment
#'
#' @param alignment alignment in DNAbin, DNAStringSet, list, and matrix formats
#'
#' @param min.n.seq minimum number of sequences needed to keep external colums
#'
#' @param codon.trim trim in triplets to avoid harming reading frame
#'
#' @return returns DNAStringSet of externally trimmed alignment
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#' @export


trimExternal = function(alignment = NULL,
                        min.n.seq = 4,
                        codon.trim = T){

  options(stringsAsFactors = FALSE)

  #Converts DNAStringSet to something usable
  #alignment<-m.align
  #min.n.seq<-50
  new.align<-strsplit(as.character(alignment), "")
  mat.align<-lapply(new.align, tolower)
  x<-as.matrix(ape::as.DNAbin(mat.align))

  if (!inherits(x, "DNAbin")) {
    stop("'x' is not of class 'DNAbin'")
  }
  if (!is.matrix(x)) {
    stop("'x' must be a matrix")
  }

  replaceWithN <- function(x) {
    id <- x == as.raw(4)
    if (length(id) > 0 & any(id[c(1, length(id))])) {
      id <- which(id)
      getIndex <- function(x) {
        for (i in seq_along(id) - 1) {
          if (any(id[1:(i + 1)] != (1:(i + 1))))
            break
        }
        id <- rev(id)
        jj <- head(id, 1)
        j <- jj - 1
        for (k in seq_along(id)[-1]) {
          if (any(id[1:k] != (jj:j)))
            break
          j <- j - 1
        }
        j <- j + 1
        id <- c(0:i, j:jj)
        id[id != 0]
      }
      id <- getIndex(id)
      x[id] <- as.raw(240)
    }
    return(x)
  }

  #Does stuff
  x <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b) {
    length(which(x %in% b))
  }
  m <- apply(x, 2, percentInformation, b)
  if (max(m) < min.n.seq) { return(0) }
  m <- range(which(m >= min.n.seq))

  #Forward Frame 2
  if (codon.trim == T){
    if ((m[1]-1) %% 3 == 0){ m[1]<-m[1] }
    if ((m[1]-1) %% 3 == 1){ m[1]<-m[1]+2 }
    if ((m[1]-1) %% 3 == 2){ m[1]<-m[1]+1 }
  }

  m <- seq(from = m[1], to = m[2])
  x2 <- as.matrix(x[, m])
  #Converts back
  save.names<-rownames(x2)

  #Removes N end gaps
  x3<-as.list(data.frame(t(as.character(x2))))
  for (y in 1:length(x3)){
    #Starts from the beginning and end to fill in end gaps
    for (q in 1:length(x3[[y]])){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
    for (q in length(x3[[y]]):1){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
  }#end x loop
  #Saves final stuff
  temp.align<-lapply(x3, FUN = function(x) paste(x, collapse = ""))
  align.out<-Biostrings::DNAStringSet(unlist(temp.align))
  names(align.out)<-save.names
  return(align.out)
}


