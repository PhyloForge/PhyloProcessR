source("workflow-X5_configuration-file.R")
if (isTRUE(get0("install.latest.github", ifnotfound = FALSE))) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    stop("Install the remotes package to use install.latest.github = TRUE.")
  }
  remotes::install_github("PhyloForge/PhyloProcessR", upgrade = "never",
                          dependencies = FALSE)
}
library(PhyloProcessR)
setwd(working.directory)

dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
source.files = c(
  if (dir.exists(alignment.directory)) {
    list.files(alignment.directory, full.names = TRUE)
  } else character(),
  if (dir.exists(paralog.directory)) {
    list.files(paralog.directory, full.names = TRUE)
  } else character(),
  if (dir.exists(primary.directory)) {
    list.files(primary.directory, full.names = TRUE)
  } else character(),
  if (dir.exists(candidate.directory)) {
    list.files(candidate.directory, pattern = "_target-candidates\\.csv$",
               full.names = TRUE)
  } else character(),
  target.file, reference.file, reference.table, feature.gene.names,
  verified.copy.mapping
)
source.files = sort(unique(source.files[
  !is.na(source.files) & nzchar(source.files) & file.exists(source.files)
]))
source.manifest = data.frame(
  Path = normalizePath(source.files),
  MD5 = unname(tools::md5sum(source.files)), stringsAsFactors = FALSE
)
run.settings = list(
  target.names = sort(target.names), alignment.format = alignment.format,
  alignment.algorithm = alignment.algorithm, run.TrimAl = run.TrimAl,
  min.external.percent = min.external.percent,
  min.column.gap.percent = min.column.gap.percent,
  min.coverage.percent = min.coverage.percent,
  min.coverage.bp = min.coverage.bp,
  min.alignment.length = min.alignment.length,
  min.taxa.alignment = min.taxa.alignment,
  max.alignment.gap.percent = max.alignment.gap.percent,
  tree.model = tree.model, bootstrap.replicates = bootstrap.replicates,
  tree.seed = tree.seed, min.branch.support = min.branch.support,
  min.shared.samples = min.shared.samples,
  min.shared.sample.fraction = min.shared.sample.fraction,
  min.split.branch.length = min.split.branch.length,
  min.split.branch.ratio = min.split.branch.ratio,
  review.action = review.action, source.manifest = source.manifest
)
settings.file = file.path(output.directory, "run-settings.rds")
if (file.exists(settings.file) && !overwrite &&
    !identical(readRDS(settings.file), run.settings)) {
  stop("X5 settings changed. Use overwrite = TRUE or a new output directory.")
}
saveRDS(run.settings, settings.file)

##################################################################################################
## Step 1: Collect base alignments and candidate copies
##################################################################################################

if (collect.copies == TRUE) {
  collectParalogCopies(
    alignment.directory = alignment.directory,
    alignment.format = alignment.format,
    paralog.directory = paralog.directory,
    primary.directory = primary.directory,
    candidate.directory = candidate.directory,
    output.directory = output.directory,
    target.names = target.names,
    reference.file = reference.file,
    reference.table = reference.table,
    overwrite = overwrite
  )
}

##################################################################################################
## Step 2: Expand alignments with candidate copies
##################################################################################################

if (align.copies == TRUE) {
  alignParalogCopies(
    collected.directory = file.path(output.directory, "1_collected"),
    target.file = target.file,
    output.directory = output.directory,
    alignment.algorithm = alignment.algorithm,
    threads = threads,
    mafft.path = mafft.path,
    quiet = quiet,
    overwrite = overwrite
  )
}

##################################################################################################
## Step 3: Prepare copy-aware alignments for tree inference
##################################################################################################

if (trim.copies == TRUE) {
  trimParalogCopies(
    expanded.directory = file.path(output.directory, "2_expanded"),
    collected.directory = file.path(output.directory, "1_collected"),
    output.directory = output.directory,
    run.TrimAl = run.TrimAl,
    trimAl.path = trimAl.path,
    min.external.percent = min.external.percent,
    min.column.gap.percent = min.column.gap.percent,
    min.coverage.percent = min.coverage.percent,
    min.coverage.bp = min.coverage.bp,
    min.alignment.length = min.alignment.length,
    min.taxa.alignment = min.taxa.alignment,
    max.alignment.gap.percent = max.alignment.gap.percent,
    quiet = quiet,
    overwrite = overwrite
  )
}

##################################################################################################
## Step 4: Infer gene trees
##################################################################################################

if (infer.trees == TRUE) {
  inferParalogTrees(
    trimmed.directory = file.path(output.directory, "3_tree-alignments"),
    collected.directory = file.path(output.directory, "1_collected"),
    output.directory = output.directory,
    tree.model = tree.model,
    bootstrap.replicates = bootstrap.replicates,
    tree.seed = tree.seed,
    threads = threads,
    memory = memory,
    iqtree.path = iqtree.path,
    iqtree.executable = iqtree.executable,
    quiet = quiet,
    overwrite = overwrite
  )
}

##################################################################################################
## Step 5: Assess splits and export accepted markers
##################################################################################################

if (separate.copies == TRUE) {
  separateParalogCopies(
    collected.directory = file.path(output.directory, "1_collected"),
    expanded.directory = file.path(output.directory, "2_expanded"),
    trimmed.directory = file.path(output.directory, "3_tree-alignments"),
    tree.directory = file.path(output.directory, "4_trees"),
    output.directory = output.directory,
    min.branch.support = min.branch.support,
    min.shared.samples = min.shared.samples,
    min.shared.sample.fraction = min.shared.sample.fraction,
    min.split.branch.length = min.split.branch.length,
    min.split.branch.ratio = min.split.branch.ratio,
    review.action = review.action,
    run.TrimAl = run.TrimAl,
    trimAl.path = trimAl.path,
    min.external.percent = min.external.percent,
    min.column.gap.percent = min.column.gap.percent,
    min.coverage.percent = min.coverage.percent,
    min.coverage.bp = min.coverage.bp,
    min.alignment.length = min.alignment.length,
    min.taxa.alignment = min.taxa.alignment,
    max.alignment.gap.percent = max.alignment.gap.percent,
    quiet = quiet,
    overwrite = overwrite
  )
}

##################################################################################################
## Step 6: Build optional gene and unlinked datasets
##################################################################################################

if (build.downstream.datasets == TRUE) {
  if (is.null(feature.gene.names)) {
    stop("feature.gene.names is required for downstream X5 datasets.")
  }
  marker.table = data.table::fread(
    file.path(output.directory, "tables", "output-markers.tsv")
  )
  gene.table = data.table::fread(feature.gene.names)
  if (!all(c("marker", "gene") %in% names(gene.table))) {
    stop("feature.gene.names must contain marker and gene columns.")
  }
  unsplit.metadata = merge(
    marker.table[Copy_label == "" & Unlinked_eligible == TRUE,
                 .(marker = Output_marker, Source_target)],
    gene.table[, .(Source_target = marker, gene)],
    by = "Source_target", all = FALSE
  )[, .(marker, gene)]
  gene.metadata = unsplit.metadata
  unlinked.metadata = unsplit.metadata

  if (!is.null(verified.copy.mapping)) {
    verified = data.table::fread(verified.copy.mapping)
    required = c("Output_marker", "Gene_copy_id", "Unlinked_eligible")
    if (!all(required %in% names(verified))) {
      stop("verified.copy.mapping must contain Output_marker, Gene_copy_id, and Unlinked_eligible.")
    }
    if (any(!verified$Output_marker %in% marker.table$Output_marker)) {
      stop("verified.copy.mapping contains an unknown output marker.")
    }
    if (anyDuplicated(verified$Output_marker) ||
        any(is.na(verified$Gene_copy_id) | verified$Gene_copy_id == "")) {
      stop("verified.copy.mapping needs one nonempty Gene_copy_id per output marker.")
    }
    logical.values = toupper(as.character(verified$Unlinked_eligible))
    if (any(is.na(logical.values) |
            !logical.values %in% c("TRUE", "FALSE"))) {
      stop("verified.copy.mapping Unlinked_eligible values must be TRUE or FALSE.")
    }
    verified$Unlinked_eligible = logical.values == "TRUE"
    approved = verified[, .(marker = Output_marker, gene = Gene_copy_id)]
    gene.metadata = unique(rbind(gene.metadata, approved))
    approved.unlinked = verified[Unlinked_eligible == TRUE,
                                 .(marker = Output_marker,
                                   gene = Gene_copy_id)]
    unlinked.metadata = unique(rbind(unlinked.metadata, approved.unlinked))
  }
  if (anyDuplicated(gene.metadata$marker) ||
      anyDuplicated(unlinked.metadata$marker)) {
    stop("Each output marker can map to only one downstream gene.")
  }
  gene.metadata.file = file.path(output.directory, "tables",
                                 "downstream-gene-metadata.tsv")
  unlinked.metadata.file = file.path(output.directory, "tables",
                                     "downstream-unlinked-metadata.tsv")
  data.table::fwrite(gene.metadata, gene.metadata.file, sep = "\t")
  data.table::fwrite(unlinked.metadata, unlinked.metadata.file, sep = "\t")
  marker.table$Downstream_gene = gene.metadata$gene[
    match(marker.table$Output_marker, gene.metadata$marker)
  ]
  marker.table$Downstream_gene[is.na(marker.table$Downstream_gene)] = ""
  data.table::fwrite(
    marker.table,
    file.path(output.directory, "tables", "output-markers.tsv"), sep = "\t"
  )
  omissions = marker.table[
    !Output_marker %in% unlinked.metadata$marker,
    .(Output_marker, Source_target, Copy_label,
      Reason = "not_approved_for_unlinked_dataset")
  ]
  data.table::fwrite(
    omissions,
    file.path(output.directory, "tables", "downstream-omissions.tsv"),
    sep = "\t"
  )
  concatenateGenes(
    alignment.folder = file.path(output.directory, "trimmed_all-markers"),
    output.folder = file.path(output.directory, "trimmed_genes"),
    feature.gene.names = gene.metadata.file,
    input.format = "phylip", output.format = "phylip",
    minimum.exons = minimum.exons, remove.reverse = FALSE,
    overwrite = overwrite, threads = threads, memory = memory
  )
  unlinked.gene.directory = file.path(output.directory,
                                      "trimmed_unlinked_genes")
  concatenateGenes(
    alignment.folder = file.path(output.directory, "trimmed_all-markers"),
    output.folder = unlinked.gene.directory,
    feature.gene.names = unlinked.metadata.file,
    input.format = "phylip", output.format = "phylip",
    minimum.exons = minimum.exons, remove.reverse = FALSE,
    overwrite = overwrite, threads = threads, memory = memory
  )
  gatherUnlinked(
    gene.alignment.directory = unlinked.gene.directory,
    exon.alignment.directory = file.path(output.directory, "trimmed_all-markers"),
    output.directory = file.path(output.directory, "trimmed_all-unlinked"),
    feature.gene.names = unlinked.metadata.file, overwrite = overwrite
  )
}

### End workflow
