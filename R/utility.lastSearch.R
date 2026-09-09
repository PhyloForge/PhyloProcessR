# Internal helpers for the LAST searches. These functions are not exported.
#
# LAST replaces blastn for the probe matching steps. blastn cannot match a contig
# to a probe past about 25 percent divergence, and sequence capture works to
# about 35 percent. On a test of 20 contigs at 35 percent divergence, blastn
# matched 7 and LAST matched all 20, with no wrong pairings. Numbers in
# HANDOFF.md.


# Builds a LAST database from a reference file.
.lastBuildDB = function(reference.file = NULL,
                        db.prefix = NULL,
                        lastdb.command = NULL,
                        threads = 1,
                        quiet = TRUE) {

  .runCommand(paste0(lastdb.command, " -P ", threads, " ",
                     shQuote(db.prefix), " ", shQuote(reference.file)),
              quiet = quiet, task = "LAST database")

  return(invisible(NULL))
}#end .lastBuildDB


# Runs a LAST search and writes the hits to out.file.
#
# The BlastTab+ format gives the same first 14 columns as the BLAST format of
# this package, so a caller reads the file with its own headers. The last column
# is the raw score in place of the gap count, and no filter uses that column.
#
# last-train is not used. It learns one substitution rate for the whole input,
# and most contigs are near-identical to their target. The trained matrix then
# rejects the few divergent contigs, which are the ones LAST must find. Without
# training LAST finds both.
.lastSearch = function(query.file = NULL,
                       db.prefix = NULL,
                       out.file = NULL,
                       lastal.command = NULL,
                       threads = 1,
                       quiet = TRUE) {

  raw.file = paste0(out.file, ".last-raw-", Sys.getpid())
  on.exit(unlink(raw.file), add = TRUE)

  # Keep the LAST status separate from the header filtering status. A failed
  # executable must not be mistaken for a successful search with no hits.
  .runCommand(paste0(lastal.command, " -P ", threads, " -f BlastTab+ ",
                     shQuote(db.prefix), " ", shQuote(query.file),
                     " > ", shQuote(raw.file)),
              quiet = quiet, task = "LAST search", keep.stdout = TRUE)

  grep.status = suppressWarnings(system(paste0("grep -v '^#' ", shQuote(raw.file),
                                                " > ", shQuote(out.file))))
  if (grep.status == 1) {
    file.create(out.file)
  } else if (grep.status != 0) {
    stop("The LAST result filtering step failed with exit status ", grep.status, ".")
  }

  return(invisible(NULL))
}#end .lastSearch


# Adds a trailing slash to a program directory, or returns an empty string. This
# matches the path handling of the other functions in the package.
.programPrefix = function(program.path = NULL) {

  if (is.null(program.path) == TRUE || nchar(program.path) == 0) return("")

  return(paste0(sub("/+$", "", program.path), "/"))
}#end .programPrefix
