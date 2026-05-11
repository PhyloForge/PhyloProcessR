#' @title assemblyStats
#'
#' @description Computes basic summary statistics for a directory of genome or
#'   contig assembly FASTA files. For each sample, reports total nucleotides,
#'   total megabases, total contig count, and mean, median, maximum, and minimum
#'   contig lengths. Results are written to a tab-delimited summary file in the
#'   output directory.
#'
#' @param assembly.directory path to a directory containing assembly FASTA files
#'   (.fa or .fasta), one file per sample.
#'
#' @param output.directory path to the directory where the summary file
#'   \code{sample-assembly-summary.txt} will be written. Created if it does not
#'   exist. Default: \code{"assembly-stats"}.
#'
#' @param overwrite logical; if \code{TRUE} the output directory is deleted and
#'   recreated before running. Default: \code{FALSE}.
#'
#' @return Invisibly returns nothing. Writes
#'   \code{<output.directory>/sample-assembly-summary.txt}, a tab-delimited
#'   table with one row per sample and columns for total nucleotides, total
#'   megabases, total contigs, and contig length statistics.
#'
#' @export

assemblyStats = function(assembly.directory = NULL,
                         output.directory = "assembly-stats",
                         overwrite = FALSE) {


  # assembly.directory = "/Volumes/LaCie/VenomCap/data-analysis/draft-assemblies"
  # output.directory = "/Volumes/LaCie/VenomCap/data-analysis/assembly-stats"
  #
  # samtools.path = "/Users/chutter/Bioinformatics/anaconda3/envs/mitocap/bin/"
  # bwa.path = "/Users/chutter/Bioinformatics/anaconda3/envs/mitocap/bin/"
  # picard.path = "/Users/chutter/Bioinformatics/anaconda3/envs/mitocap/bin/"
  #
  # quiet = TRUE
  # overwrite = FALSE
  # threads = 6
  # memory = 6

  if (is.null(assembly.directory) == T){ stop("A directory of genome(s) is needed.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  assembly.files = list.files(assembly.directory)
  sample.names = gsub(".fa|.fasta", "", assembly.files)

  #Stats table prepare
  header.data = c("sample", "total_nucleotides", "total_megabases", "total_contigs",
                  "mean_contig_length","median_contig_length", "max_contig_length", "min_contig_length")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(assembly.files), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, sample:=as.character(sample)]

  for (i in 1:length(assembly.files)){

    assembly.data = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", assembly.files[i]), format = "fasta")

    data.table::set(collect.data, i = as.integer(i), j = match("sample", header.data), value = gsub(".fa|.fasta", "", assembly.files[i]))
    data.table::set(collect.data, i = as.integer(i), j = match("total_nucleotides", header.data), value = sum(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("total_megabases", header.data), value = sum(Biostrings::width(assembly.data))/1000000)
    data.table::set(collect.data, i = as.integer(i), j = match("total_contigs", header.data), value = length(Biostrings::width(assembly.data)) )
    data.table::set(collect.data, i = as.integer(i), j = match("mean_contig_length", header.data), value = mean(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("median_contig_length", header.data), value = median(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("max_contig_length", header.data), value = max(Biostrings::width(assembly.data)))
    data.table::set(collect.data, i = as.integer(i), j = match("min_contig_length", header.data), value = min(Biostrings::width(assembly.data)))

  }# end i loop

  write.table(collect.data, file = paste0(output.directory, "/sample-assembly-summary.txt"), sep = "\t", row.names = F)

}#end funtion


#END SCRIPT

