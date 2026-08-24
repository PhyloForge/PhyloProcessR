#' @title alignMACSE
#'
#' @description Refines and aligns exon sequences using MACSE to ensure proper reading frames and codon alignment.
#'
#' @param alignment.folder character string; path to the folder containing the alignments.
#' @param output.folder character string; path to the folder to save the MACSE-refined alignments.
#' @param macse.path character string; system path to the directory containing the macse executable. If NULL, searches the system PATH.
#' @param genetic.code integer; the genetic code table to use (default: 1 for standard nuclear).
#' @param threads integer; number of threads to use.
#' @param memory integer; memory allocated (in GB).
#' @param overwrite logical; if TRUE, overwrite existing output files.
#' @param quiet logical; if TRUE, suppress output messages.
#'
#' @return A folder of MACSE-refined codon alignments.
#'
#' @export

alignMACSE = function(alignment.folder = NULL,
                      output.folder = NULL,
                      macse.path = NULL,
                      genetic.code = 1,
                      threads = 1,
                      memory = 4,
                      overwrite = FALSE,
                      quiet = TRUE) {

  if (is.null(alignment.folder)) { stop("An input alignment folder must be provided.") }
  if (is.null(output.folder)) { stop("An output folder must be provided.") }

  # Creates output directory if it doesn't exist
  if (dir.exists(output.folder) == FALSE) {
    dir.create(output.folder, recursive = TRUE)
  }

  if (is.null(macse.path) == FALSE){
    b.string = unlist(strsplit(macse.path, ""))
    if (b.string[length(b.string)] != "/") {
      macse.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { macse.path = "" }

  # Gets all alignment files
  align.files = list.files(alignment.folder, full.names = TRUE, pattern = "\\.phy$|\\.fa$|\\.fasta$")
  if (length(align.files) == 0) { stop("No alignments found in the input folder.") }

  # Sets up foreach loop
  require(foreach)
  require(doParallel)
  
  cl = makeCluster(threads)
  registerDoParallel(cl)

  cat(paste0("Refining ", length(align.files), " alignments using MACSE...\n"))

  foreach(i = 1:length(align.files), .packages = c("Biostrings", "ape")) %dopar% {
    
    file.name = basename(align.files[i])
    file.base = strsplit(file.name, "\\.")[[1]][1]
    
    out.file = paste0(output.folder, "/", file.base, ".fa")
    
    if (file.exists(out.file) && overwrite == FALSE) {
      return(NULL)
    }

    # MACSE requires java. The bioconda macse wrapper can be called directly
    # -prog refineAlignment requires -align. -prog alignSequences requires -seq.
    # To be safe against gappy untrimmed sequences, alignSequences is typically best,
    # but since these are already aligned exons, refineAlignment is also good. 
    # Let's use alignSequences to do a fresh robust codon alignment, or refineAlignment.
    # We will use alignSequences since it cleans things up beautifully.
    macse_cmd = paste0(
      macse.path, "macse -prog alignSequences ",
      "-seq ", align.files[i], " ",
      "-gc_def ", genetic.code, " ",
      "-out_NT ", out.file, " ",
      "-out_AA /dev/null"
    )
    
    if (quiet == TRUE) {
      macse_cmd = paste0(macse_cmd, " > /dev/null 2>&1")
    }
    
    system(macse_cmd)
  }

  stopCluster(cl)
  cat("MACSE refinement complete.\n")
}
