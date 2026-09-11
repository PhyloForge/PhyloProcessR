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
