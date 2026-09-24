.alignmentExtensions = function(format = NULL) {
  if (is.null(format)) return(c("phy", "phylip", "fa", "fas", "fasta", "nex", "nexus"))

  format = match.arg(format, c("phylip", "fasta", "nexus"))
  switch(format,
         phylip = c("phy", "phylip"),
         fasta = c("fa", "fas", "fasta"),
         nexus = c("nex", "nexus"))
}

.alignmentFiles = function(directory, format = NULL, full.names = FALSE) {
  if (!dir.exists(directory)) return(character())
  extensions = paste(.alignmentExtensions(format), collapse = "|")
  list.files(directory,
             pattern = paste0("\\.(", extensions, ")$"),
             full.names = full.names,
             recursive = FALSE,
             ignore.case = TRUE)
}

.alignmentId = function(path) {
  extensions = paste(.alignmentExtensions(), collapse = "|")
  sub(paste0("\\.(", extensions, ")$"), "", basename(path),
      ignore.case = TRUE)
}

.copyAlignment = function(source, destination, overwrite = FALSE) {
  if (file.exists(destination) && !overwrite) return(FALSE)
  copied = file.copy(source, destination, overwrite = overwrite)
  if (!isTRUE(copied)) {
    stop("Could not copy alignment from ", source, " to ", destination, ".")
  }
  TRUE
}

.writeAlignmentFailureLogs = function(results, id.field, log.directory, suffix) {
  failures = vapply(results, function(x) identical(x$status, "error"), logical(1))
  if (!any(failures)) return(invisible(0L))

  dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
  for (result in results[failures]) {
    alignment.id = gsub("[/\\\\]", "_", as.character(result[[id.field]]))
    log.file = file.path(log.directory, paste0(alignment.id, suffix))
    writeLines(result$message, log.file)
  }
  invisible(sum(failures))
}

.copyAlignmentsWithLogs = function(sources, output.directory, overwrite,
                                   log.directory, suffix) {
  copied = logical(length(sources))
  for (i in seq_along(sources)) {
    destination = file.path(output.directory, basename(sources[i]))
    if (!overwrite && file.exists(destination)) {
      copied[i] = TRUE
      next
    }
    copied[i] = tryCatch(
      .copyAlignment(sources[i], destination, overwrite = overwrite),
      error = function(e) {
        dir.create(log.directory, recursive = TRUE, showWarnings = FALSE)
        alignment.id = gsub("[/\\\\]", "_", .alignmentId(sources[i]))
        writeLines(conditionMessage(e),
                   file.path(log.directory, paste0(alignment.id, suffix)))
        FALSE
      }
    )
  }
  invisible(copied)
}

.writePhylipAtomic = function(alignment, destination, interleave = FALSE,
                              strict = FALSE) {
  temporary = tempfile(paste0(".", basename(destination), "-"),
                       tmpdir = dirname(destination), fileext = ".phy")
  on.exit(unlink(temporary), add = TRUE)
  PhyloProcessR::writePhylip(alignment, file = temporary,
                             interleave = interleave, strict = strict)
  if (!file.exists(temporary) || file.info(temporary)$size == 0) {
    stop("A complete PHYLIP alignment was not written: ", destination)
  }
  Biostrings::readDNAMultipleAlignment(temporary, format = "phylip")
  if (file.exists(destination)) unlink(destination)
  if (!file.rename(temporary, destination)) {
    stop("Could not publish the PHYLIP alignment: ", destination)
  }
  invisible(destination)
}

.readGeneMetadata = function(path) {
  if (is.null(path) || !file.exists(path)) {
    stop("The gene metadata file was not found: ", path)
  }

  metadata = data.table::fread(path, header = TRUE)
  if (!"marker" %in% names(metadata) && "Marker" %in% names(metadata)) {
    data.table::setnames(metadata, "Marker", "marker")
  }
  if (!"gene" %in% names(metadata) && "Gene" %in% names(metadata)) {
    data.table::setnames(metadata, "Gene", "gene")
  }
  missing.columns = setdiff(c("marker", "gene"), names(metadata))
  if (length(missing.columns) > 0) {
    stop("Gene metadata is missing required column(s): ",
         paste(missing.columns, collapse = ", "), ".")
  }

  assigned = metadata[!is.na(gene) & nzchar(as.character(gene))]
  conflicts = assigned[, .(gene.count = data.table::uniqueN(gene)), by = marker]
  conflicts = conflicts[gene.count > 1]
  if (nrow(conflicts) > 0) {
    examples = paste(utils::head(conflicts$marker, 5), collapse = ", ")
    stop("Markers are assigned to more than one gene: ", examples, ".")
  }
  metadata
}
