#' @title alignMACSE
#'
#' @description Refines and aligns exon sequences using MACSE to ensure proper reading frames and codon alignment.
#'
#' @param alignment.folder character string; path to the folder containing the alignments.
#' @param output.folder character string; path to the folder to save the MACSE-refined alignments.
#' @param alignment.format character string; the format of the input alignments (e.g. "phylip", "fasta").
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

  if (is.null(macse.path) == FALSE){
    b.string = unlist(strsplit(macse.path, ""))
    if (b.string[length(b.string)] != "/") {
      macse.path = paste0(append(b.string, "/"), collapse = "")
    }
  } else { macse.path = "" }

  # Gets all alignment files
  if (alignment.format == "phylip") {
    align.files = list.files(alignment.folder, full.names = TRUE, pattern = "\\.phy$")
  } else {
    align.files = list.files(alignment.folder, full.names = TRUE, pattern = "\\.fa$|\\.fasta$")
  }

  if (length(align.files) == 0) { stop("No alignments found in the input folder.") }

  # Sets up foreach loop
  require(foreach)
  require(doParallel)
  
  cl = makeCluster(threads)
  registerDoParallel(cl)

  cat(paste0("Refining ", length(align.files), " alignments using MACSE...\n"))

  foreach(i = 1:length(align.files), .packages = c("Biostrings", "ape", "seqinr")) %dopar% {
    
    file.name = basename(align.files[i])
    file.base = strsplit(file.name, "\\.")[[1]][1]
    
    out.file = paste0(output.folder, "/", file.base, ".fa")
    out.aa.file = paste0(output.folder, "/temp_", file.base, "_AA.fa")
    
    if (file.exists(out.file) && overwrite == FALSE) {
      return(NULL)
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
      macse.path, "macse -prog alignSequences ",
      "-seq ", macse_input, " ",
      "-gc_def ", genetic.code, " ",
      "-out_NT ", out.file, " ",
      "-out_AA ", out.aa.file
    )
    
    log_file = paste0(output.folder, "/temp_macse_error_log.txt")
    if (quiet == TRUE) {
      macse_cmd = paste0(macse_cmd, " > ", log_file, " 2>&1")
    }
    
    system(macse_cmd)
    
    # Format conversion if needed
    if (output.format == "phylip" && file.exists(out.file)) {
      align_macse = ape::read.FASTA(out.file, type = "DNA")
      align_mat = as.matrix(align_macse)
      rownames(align_mat) = labels(align_macse)
      
      out.phy = paste0(output.folder, "/", file.base, ".phy")
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
  }

  stopCluster(cl)
  cat("MACSE refinement complete.\n")
}
