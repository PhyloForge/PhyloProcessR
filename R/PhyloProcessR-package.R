#' PhyloProcessR: A Modular Toolkit for Targeted Sequence-Capture Phylogenomics
#'
#' PhyloProcessR provides independently callable R functions and complete
#' workflow scripts for processing targeted sequence-capture data from raw
#' reads to curated phylogenomic datasets. Functions can be composed into
#' reproducible, project-specific pipelines. See the project README and
#' workflow configuration files for the external tools required by each stage.
#'
#' @keywords internal
#' @import data.table
#' @import ggplot2
#' @importFrom ape as.DNAbin getMRCA is.monophyletic read.nexus read.tree root
#' @importFrom Biostrings DNAStringSet reverseComplement subseq
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach %dopar% foreach
#' @importFrom grDevices dev.off pdf
#' @importFrom grid unit
#' @importFrom parallel makeCluster stopCluster
#' @importFrom Rsamtools FaFile scanFa
#' @importFrom stats aggregate heatmap median sd setNames
#' @importFrom stringr str_locate_all
#' @importFrom utils capture.output download.file head read.csv read.table tail
#'   write.csv write.table
"_PACKAGE"
