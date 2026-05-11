#' @title writePhylip
#'
#' @description Writes a multiple sequence alignment to phylip format. The
#'   alignment is expected to be a matrix-like object with row names used as
#'   taxon labels. Supports both sequential and interleaved output, and optionally
#'   truncates taxon names to 10 characters in strict phylip mode.
#'
#' @param alignment a matrix or matrix-convertible object (e.g. a DNAbin or
#'   character matrix) where rows are taxa and columns are alignment positions.
#'   Row names are used as taxon labels.
#'
#' @param file path to the output phylip file; if an empty string (""), the
#'   phylip text is printed to the console instead.
#'
#' @param interleave logical or integer; if FALSE (default) sequences are
#'   written sequentially. If an integer, sequences are interleaved with that
#'   number of characters per block.
#'
#' @param strict logical; if TRUE taxon names are truncated to 10 characters and
#'   padded with asterisks to produce strict phylip format. A warning is issued
#'   if truncation creates duplicate names. If FALSE names are padded with
#'   spaces to the length of the longest name plus 3.
#'
#' @return invisibly; writes the alignment to file in phylip format, or prints
#'   to the console if file is "".
#'
#' @export

##### saves alignment as a phylip file
writePhylip = function(alignment = NULL,
                       file = NULL,
                       interleave = FALSE,
                       strict = FALSE) {

  x = as.matrix(alignment)
  #file = "test"

  #Parameter checks
  if(is.null(alignment) == TRUE){ stop("Error: an input alignment is needed.") }
  if(is.null(file) == TRUE){ stop("Error: a file name is needed.") }
  if(is.null(row.names(alignment)) == TRUE){ stop("Error: no row names found on matrix.") }

  #converts alignment
  # new.align = alignmentConversion(input.alignment = alignment, end.format = "matrix")
  #
  # ntax = nrow(new.align)
  # nchar = ncol(new.align)
  # phy.header = paste(ntax, nchar)

  # #Now have to reformat somehow
  # #Prep sample names to all have the same length padded with spaces
  # name.length = max(nchar(sample.names)) + 5
  # for (x in 1:length(sample.names)){
  #   temp.name = concat.data$Sample[x]
  #   space.pad = name.length - nchar(temp.name)
  #   space.add = paste(rep(" ", space.pad), collapse = "")
  #   new.name = paste0(temp.name, space.add, collapse = "")
  #   concat.data$Sample[x] = new.name
  # }#end x

  str2cha <- function(x) { unlist(strsplit(x, "")) }
  datatype <- ifelse(is.numeric(x[1, 1]), "continuous", "nc")
  ntax <- nrow(x)
  nchar <- ncol(x)
  taxnames <- rownames(x)
  if (strict) {
    taxnames <- substring(taxnames, 1, truncate)
    missing <- 10 - unlist(lapply(strsplit(taxnames, ""),
                                  length))
    for (i in seq(along = taxnames)) taxnames[i] <- paste(taxnames[i],
                                                          paste(rep("*", missing[i]), collapse = ""), sep = "")
    if (any(duplicated(taxnames)))
      cat("WARNING: Truncation of taxon names created",
          "identical strings.")
  }
  else {
    xx <- nchar(taxnames)
    diff <- max(xx) - xx + 3
    for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], paste(rep(" ",
                                                                  diff[i]), collapse = ""), sep = "")
  }#end else

  if (!interleave)
    interleave <- nchar
  nbpart <- ceiling(nchar/interleave)
  pt <- matrix(nrow = nbpart, ncol = 2)
  pt[1, ] <- c(1, interleave)
  if (nbpart > 1)
    for (i in 2:(dim(pt)[1])) {
      pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
      pt[nbpart, 2] <- nchar
    }#end for loop

  phy <- paste(ntax, nchar)
  for (i in seq(along = pt[, 1])) {
    sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
    if (is.null(dim(sm)))
      sm <- as.matrix(sm, ncol = 1)
    sm <- apply(sm, 1, paste, collapse = "")
    if (i == 1)
      sm <- paste(taxnames, sm)
    if (i < max(seq(along = pt[, 1])))
      sm <- c(sm, "")
    phy <- c(phy, sm)
  } #end for loop
  if (file == "") {
    cat(phy, sep = "\n")
  }#end if
  else {
    write(phy, file = file)
  }#end else

}
