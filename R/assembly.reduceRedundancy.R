#' @title reduceRedundancy
#'
#' @description Clusters and de-replicates contigs within each sample assembly
#'   using \code{cd-hit-est} at a user-specified sequence identity threshold.
#'   Each FASTA file in the input directory is processed in parallel: cd-hit-est
#'   collapses highly similar sequences, contigs are renamed sequentially, and
#'   the dereplicated assembly is written to the output directory. The cluster
#'   report (\code{.clstr}) files are removed after processing.
#'
#' @param assembly.directory path to a directory of assembly FASTA files, one
#'   per sample.
#'
#' @param output.directory path to the directory where dereplicated FASTA files
#'   will be written. Default: \code{"reduced-redundancy"}.
#'
#' @param similarity sequence identity threshold for cd-hit-est clustering
#'   (0.8-1.0). Must be >= 0.8, which is the lowest threshold cd-hit-est
#'   accepts. Default: \code{0.95}.
#'
#' @param cdhit.path path to the directory containing \code{cd-hit-est}. If
#'   \code{NULL} expected on the system PATH. Default: \code{NULL}.
#'
#' @param memory total RAM in GB to divide across parallel workers. Default:
#'   \code{1}.
#'
#' @param threads number of parallel workers. Default: \code{1}.
#'
#' @param overwrite logical; if \code{TRUE} the output directory is deleted and
#'   recreated. Default: \code{FALSE}.
#'
#' @param quiet logical; if \code{TRUE} cd-hit-est screen output is written to
#'   the per-sample log instead of the console. Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Writes one dereplicated FASTA file per
#'   sample to \code{output.directory}, with contigs renamed
#'   \code{contig_1}, \code{contig_2}, etc.
#'
#' @export

#Iteratively assembles to reference
reduceRedundancy = function(assembly.directory = NULL,
                            output.directory = "reduced-redundancy",
                            similarity = 0.95,
                            cdhit.path = NULL,
                            memory = 1,
                            threads = 1,
                            overwrite = FALSE,
                            quiet = TRUE) {

  # #Debug
  # library(PhyloProcessR)
  # setwd("/Volumes/LaCie/Mantellidae")
  # assembly.directory = "data-analysis/contigs/draft-assemblies"
  # output.directory = "data-analysis/contigs/reduced-redundancy"
  # cdhit.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # similarity <- 0.95
  # quiet = TRUE
  # overwrite = TRUE
  # threads = 5
  # memory = 30

  # Quick checks
  if (is.null(assembly.directory) == TRUE) { stop("Please provide an assembly directory.") }
  if (is.null(output.directory) == TRUE)   { stop("Please provide an output directory.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE) {
      unlink(output.directory, recursive = TRUE)
      dir.create(output.directory, recursive = TRUE)
    }
  } else {
    dir.create(output.directory, recursive = TRUE)
  }

  # Only FASTA files are samples. Hidden files and leftover .clstr reports are
  # ignored.
  fasta.pattern = "\\.fa$|\\.fas$|\\.fasta$|\\.fna$"
  file.names = list.files(assembly.directory, pattern = fasta.pattern)

  # Resume: skip samples already written to the output directory
  if (overwrite == FALSE) {
    done = list.files(output.directory, pattern = fasta.pattern)
    file.names = file.names[!file.names %in% done]
  }

  if (length(file.names) == 0) { return(invisible(NULL)) }

  if (!is.numeric(threads) || length(threads) != 1 || !is.finite(threads) || threads < 1 ||
      !is.numeric(memory) || length(memory) != 1 || !is.finite(memory) || memory <= 0) {
    stop("threads and memory must be positive finite values.")
  }
  workers = min(floor(threads), length(file.names), max(1, floor(memory)))
  cdhit.command = .toolCommand("cd-hit-est", cdhit.path)

  # Word length table from the cd-hit-est manual. cd-hit-est stops with a fatal
  # error below a threshold of 0.8.
  if (similarity >= 0.95) {
    n.val = 10
  } else if (similarity >= 0.90) {
    n.val = 8
  } else if (similarity >= 0.88) {
    n.val = 7
  } else if (similarity >= 0.85) {
    n.val = 6
  } else if (similarity >= 0.80) {
    n.val = 5
  } else {
    stop("similarity too small for cd-hit-est. Must be >= 0.8.")
  }

  # Memory is split across the workers. A value of 0 means unlimited in cd-hit,
  # so the split is never allowed to reach 0.
  mem.cl = max(1, floor(memory / workers))

  if (dir.exists("logs/sample_logs") == FALSE) { dir.create("logs/sample_logs", recursive = TRUE) }

  # Use mclapply (fork-based) instead of a SOCK cluster.
  # Fork workers die silently per-sample rather than crashing the whole session
  # when a socket connection is lost -- the primary cause of the previous crash.
  parallel::mclapply(seq_along(file.names), function(i) {

    sample.id  = file.names[i]
    out.file   = paste0(output.directory, "/", sample.id)
    # cd-hit-est writes to a temporary name. The final file only appears after the
    # contigs are renamed, so an interrupted run is not skipped on the next resume.
    tmp.file   = paste0(out.file, ".tmp")
    log.file   = paste0("logs/sample_logs/FAILURE_", sample.id, "_reduceRedundancy.txt")
    cdhit.log  = paste0("logs/sample_logs/", sample.id, "_cdhit-output.txt")

    tryCatch({

      # Check input is readable before handing to cd-hit
      in.file = paste0(assembly.directory, "/", sample.id)
      if (!file.exists(in.file) || file.info(in.file)$size == 0) {
        msg = paste0("Input FASTA missing or empty: ", in.file)
        writeLines(msg, log.file)
        warning(sample.id, ": ", msg)
        return(invisible(NULL))
      }

      # Redirect cd-hit-est output to a per-sample log so we can inspect failures.
      # Parallel workers would otherwise interleave their screen output.
      exit.code = system(paste0(
        cdhit.command, " -i ", shQuote(in.file),
        " -o ", shQuote(tmp.file), " -p 0 -T 1",
        " -n ", n.val, " -c ", similarity, " -M ", mem.cl * 1000,
        " > ", shQuote(cdhit.log), " 2>&1"
      ))

      if (exit.code != 0) {
        cdhit.msg = if (file.exists(cdhit.log)) paste(readLines(cdhit.log, warn = FALSE), collapse = "\n") else "(no output captured)"
        msg = paste0("cd-hit-est exited with code ", exit.code, ".\n\ncd-hit-est output:\n", cdhit.msg)
        writeLines(msg, log.file)
        warning(sample.id, ": cd-hit-est failed (exit code ", exit.code,
                ") -- see ", log.file)
        unlink(c(tmp.file, paste0(tmp.file, ".clstr")))
        return(invisible(NULL))
      }

      if (!file.exists(tmp.file) || file.info(tmp.file)$size == 0) {
        msg = "cd-hit-est exited 0 but produced no output file."
        writeLines(msg, log.file)
        warning(sample.id, ": ", msg)
        unlink(c(tmp.file, paste0(tmp.file, ".clstr")))
        return(invisible(NULL))
      }

      all.data = Biostrings::readDNAStringSet(file = tmp.file, format = "fasta")

      if (length(all.data) == 0) {
        msg = "cd-hit-est output FASTA contains no sequences."
        writeLines(msg, log.file)
        warning(sample.id, ": ", msg)
        unlink(c(tmp.file, paste0(tmp.file, ".clstr")))
        return(invisible(NULL))
      }

      names(all.data) = paste0("contig_", seq_along(all.data))

      final.loci = as.list(as.character(all.data))
      publish.file = paste0(out.file, ".publish-", Sys.getpid())
      PhyloProcessR::writeFasta(
        sequences = final.loci, names = names(final.loci),
        publish.file,
        nbchar = 1000000, as.string = TRUE, open = "w"
      )
      if (!file.exists(publish.file) || file.size(publish.file) == 0 ||
          file.rename(publish.file, out.file) == FALSE) {
        unlink(publish.file)
        stop("Could not publish reduced contigs for ", sample.id, ".")
      }

      unlink(c(tmp.file, paste0(tmp.file, ".clstr")))

      # Clean up the output log on success -- only keep it for failures
      if (file.exists(cdhit.log)) {
        if (quiet == FALSE) { cat(readLines(cdhit.log, warn = FALSE), sep = "\n") }
        file.remove(cdhit.log)
      }

    }, error = function(e) {
      msg = paste0("Unexpected R error: ", conditionMessage(e), "\n\n",
                   paste(capture.output(traceback()), collapse = "\n"))
      writeLines(msg, log.file)
      warning(sample.id, ": unexpected error -- see ", log.file)
    })

  }, mc.cores = workers)

  invisible(NULL)

}#end function




