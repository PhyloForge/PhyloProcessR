#' @title alignMACSE
#'
#' @description Refines and aligns exon sequences using MACSE to ensure proper reading frames and codon alignment.
#'
#' @param alignment.folder character string; path to the folder containing the alignments.
#' @param output.folder character string; path to the folder to save the MACSE-refined alignments.
#' @param alignment.format character string; the format of the input alignments (e.g. "phylip", "fasta").
#' @param output.format character string; output alignment format. Currently
#'   used to name the intended output format; PHYLIP output is supported.
#' @param macse.path character string; system path to the directory containing the macse executable. If NULL, searches the system PATH.
#' @param genetic.code integer; the genetic code table to use (default: 1 for standard nuclear).
#' @param threads integer; number of threads to use.
#' @param memory integer; reserved for compatibility. This value does not set the JVM memory limit.
#' @param overwrite logical; if TRUE, overwrite existing output files.
#' @param quiet logical; if TRUE, suppress output messages.
#'
#' @return A folder of MACSE-refined codon alignments.
#'
#' @export

alignMACSE = function(alignment.folder = NULL,
                      output.folder = NULL,
                      alignment.format = "phylip",
                      output.format = "phylip",
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

  macse.command = if (is.null(macse.path)) "macse" else
    file.path(macse.path, "macse")

  # Gets all alignment files
  if (alignment.format == "phylip") {
    align.files = list.files(alignment.folder, full.names = TRUE, pattern = "\\.phy$")
  } else {
    align.files = list.files(alignment.folder, full.names = TRUE, pattern = "\\.fa$|\\.fasta$")
  }

  if (length(align.files) == 0) { stop("No alignments found in the input folder.") }

  # Sets up foreach loop
  
  cl = makeCluster(threads)
  registerDoParallel(cl)
  on.exit({
    stopCluster(cl)
    foreach::registerDoSEQ()
  }, add = TRUE)

  cat(paste0("Refining ", length(align.files), " alignments using MACSE...\n"))

  results = foreach(i = seq_along(align.files),
                    .packages = c("Biostrings", "ape", "seqinr")) %dopar% {
    
    file.name = basename(align.files[i])
    file.base = .alignmentId(file.name)
    
    out.file = paste0(output.folder, "/", file.base, ".fa")
    out.aa.file = paste0(output.folder, "/temp_", file.base, "_AA.fa")
    
    final.file = if (output.format == "phylip") {
      file.path(output.folder, paste0(file.base, ".phy"))
    } else {
      out.file
    }
    if (file.exists(final.file) && file.info(final.file)$size > 0 && !overwrite) {
      return(list(status = "skipped", locus = file.base))
    }

    # MACSE requires FASTA format. If input is phylip, we need to convert to a temp fasta file.
    macse_input = align.files[i]
    
    if (alignment.format == "phylip") {
      align = ape::read.dna(align.files[i], format = "sequential")
      temp.align = as.character(as.list(align))
      temp.align2 = lapply(temp.align, FUN = function(x) paste(x, collapse = ""))
      align.out = Biostrings::DNAStringSet(unlist(temp.align2))
      
      temp.fa = paste0(output.folder, "/temp_", file.base, ".fa")
      
      write.loci = as.list(as.character(align.out))
      seqinr::write.fasta(sequences = write.loci, names = names(write.loci),
                          file.out = temp.fa, nbchar = 1000000, as.string = TRUE)
                          
      macse_input = temp.fa
    }

    # MACSE command
    macse_cmd = paste0(
      shQuote(macse.command), " -prog alignSequences ",
      "-seq ", shQuote(macse_input), " ",
      "-gc_def ", genetic.code, " ",
      "-out_NT ", shQuote(out.file), " ",
      "-out_AA ", shQuote(out.aa.file)
    )
    
    log.directory = file.path(output.folder, "logs")
    dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
    log_file = file.path(log.directory, paste0(file.base, "_macse.log"))
    status = system(paste0(macse_cmd, " > ", shQuote(log_file), " 2>&1"))
    if (status != 0) {
      return(list(status = "error", locus = file.base,
                  message = paste0("MACSE exited with status ", status, ".")))
    }
    if (!file.exists(out.file) || file.info(out.file)$size == 0) {
      return(list(status = "error", locus = file.base,
                  message = "MACSE did not create a nucleotide alignment."))
    }
    
    # Format conversion if needed
    if (output.format == "phylip" && file.exists(out.file)) {
      align_macse = ape::read.FASTA(out.file, type = "DNA")
      align_mat = as.matrix(align_macse)
      rownames(align_mat) = labels(align_macse)
      
      out.phy = final.file
      PhyloProcessR::writePhylip(align_mat, file = out.phy)
      
      # Delete the fasta output from MACSE
      file.remove(out.file)
    }
    
    # Cleanup temp files
    if (alignment.format == "phylip" && file.exists(temp.fa)) {
      file.remove(temp.fa)
    }
    if (file.exists(out.aa.file)) {
      file.remove(out.aa.file)
    }
    list(status = "success", locus = file.base)
  }

  failures = vapply(results, function(x) identical(x$status, "error"), logical(1))
  if (any(failures)) {
    details = vapply(results[failures], function(x) {
      paste0(x$locus, " (", x$message, ")")
    }, character(1))
    stop("MACSE failed for: ", paste(details, collapse = "; "))
  }
  cat("MACSE refinement complete.\n")
}
