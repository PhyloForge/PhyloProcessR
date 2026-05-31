#' @title trimTrimal
#'
#' @description Wrapper function for running TrimAl on a single alignment. The alignment is written to a temporary fasta file, TrimAl is called in automated mode (-automated1), and the trimmed alignment is read back. Sample names are restored from the original alignment to correct any truncation introduced by TrimAl. If TrimAl produces no output the original alignment is returned unchanged. Alignments with three or fewer sequences are returned unmodified. TrimAl must be installed and accessible.
#'
#' @param alignment a DNAStringSet containing the aligned sequences to trim
#'
#' @param trimal.path system path to the directory containing the trimal executable; NULL to use the system PATH
#'
#' @param method trimming method to use; currently only "auto" (automated1) is implemented
#'
#' @param quiet if TRUE, suppress TrimAl screen output
#'
#' @return a DNAStringSet of the TrimAl-trimmed alignment with original sample names restored; returns the original alignment unchanged if TrimAl produces no output
#'
#' @export

trimTrimal = function(alignment = NULL,
                      trimal.path = NULL,
                      method = "auto",
                      quiet = TRUE) {
  #Debug
   # alignment = non.align
   # quiet = FALSE
   # trimal.path = "/Users/chutter/conda/PhyloCap/bin"

  #Same adds to bbmap path
  if (is.null(trimal.path) == FALSE){
    b.string = unlist(strsplit(trimal.path, ""))
    if (b.string[length(b.string)] != "/") {
      trimal.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { trimal.path = "" }

  if (length(alignment) <= 3){ return(alignment) }

  #Finds probes that match to two or more contigs
  save.rownames = names(alignment)
  write.align = as.list(as.character(alignment))

  #Creates random name and saves it
  input.file = paste0("temp_", sample(1:1000000, 1), ".fa")
  writeFasta(sequences = write.align,
             names = names(write.align),
             file.out = input.file,
             nbchar = 1000000,
             as.string = T)

  #Runs trimal command with input file
  system(paste0(trimal.path, "trimal -in ", input.file, " -out tm-", input.file, " -automated1"),
         ignore.stdout = quiet, ignore.stderr = quiet)

  if (file.exists(paste0("tm-", input.file)) == F) {
    system(paste0("rm ", input.file))
    if (file.exists(paste0(input.file, ".fai"))) { system(paste0("rm ", input.file, ".fai")) }
    return(alignment)
  } else { system(paste0("mv tm-", input.file, " ", input.file)) }

  out.align = Rsamtools::scanFa(Rsamtools::FaFile(input.file))

  # Restore original (full-length) names from save.rownames.
  # The old approach used an unescaped regex grep with a "$" anchor, which has two failure modes:
  #   (1) Suffix collision: if name A is a suffix of name B, grep matches BOTH and the
  #       multi-element assignment silently corrupts new.names, shifting every subsequent label.
  #   (2) Regex metacharacters: dots and other special chars in taxon names cause false matches.
  # Fix: exact match first; fall back to literal endsWith() for genuinely truncated names.
  new.names = vapply(names(out.align), function(nm) {
    exact = which(save.rownames == nm)
    if (length(exact) == 1L) return(save.rownames[exact])
    suffix = which(endsWith(save.rownames, nm))
    if (length(suffix) == 1L) return(save.rownames[suffix])
    nm  # last resort: keep whatever TrimAl gave us
  }, character(1L))

  temp = names(out.align)[is.na(names(out.align)) == T]
  if (length(temp) > 0){ stop("there are NAs in the names") }
  names(out.align) = new.names
  system(paste0("rm ", input.file))
  if (file.exists(paste0(input.file, ".fai"))) { system(paste0("rm ", input.file, ".fai")) }
  return(out.align)

}#end function
