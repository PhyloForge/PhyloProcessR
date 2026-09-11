#' @title trimTrimal
#'
#' @description Wrapper function for running TrimAl on a single alignment. The alignment is written to a temporary fasta file, TrimAl is called in automated mode (-automated1), and the trimmed alignment is read back. Sample names are restored from the original alignment to correct any truncation introduced by TrimAl. Alignments with three or fewer sequences are returned unmodified. TrimAl must be installed and accessible.
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

  trimal.command = if (is.null(trimal.path)) "trimal" else
    file.path(trimal.path, "trimal")

  if (length(alignment) <= 3){ return(alignment) }

  #Finds probes that match to two or more contigs
  save.rownames = names(alignment)
  write.align = as.list(as.character(alignment))

  #Creates random name and saves it
  input.file = paste0("temp_", sample(1:1000000, 1), ".fa")
  output.file = paste0("tm-", input.file)
  log.file = paste0(input.file, ".log")
  on.exit(unlink(c(input.file, paste0(input.file, ".fai"), output.file)), add = TRUE)
  writeFasta(sequences = write.align,
             names = names(write.align),
             file.out = input.file,
             nbchar = 1000000,
             as.string = T)

  #Runs trimal command with input file
  .runCommand(paste0(shQuote(trimal.command), " -in ", shQuote(input.file),
                     " -out ", shQuote(output.file), " -automated1"),
              quiet = quiet, task = "trimAl trimming", stderr.log = log.file)
  if (!file.exists(output.file) || file.info(output.file)$size == 0) {
    stop("trimAl did not create a readable output alignment.")
  }
  if (!file.rename(output.file, input.file)) stop("Could not prepare trimAl output.")

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
  unlink(log.file)
  return(out.align)

}#end function
