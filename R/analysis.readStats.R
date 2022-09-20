#' @title readStats
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param genome.directory path to a folder of sequence alignments in phylip format.
#'
#' @param output.directory available input alignment formats: fasta or phylip
#'
#' @param threads contigs are added into existing alignment if algorithm is "add"
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite TRUE to supress mafft screen output
#'
#' @param resume TRUE to supress mafft screen output
#'
#' @return an alignment of provided sequences in DNAStringSet format. Also can save alignment as a file with save.name
#'
#' @examples
#'
#' your.tree = ape::read.tree(file = "file-path-to-tree.tre")
#' astral.data = astralPlane(astral.tree = your.tree,
#'                           outgroups = c("species_one", "species_two"),
#'                           tip.length = 1)
#'
#'
#' @export


analysis.readStats = function(read.directory = NULL,
                              output.directory = "read-stats",
                              overwrite = FALSE) {

  #Debug comment out
  # read.directory = "/Volumes/LaCie/VenomCap/read-processing/cleaned-reads"
  # output.directory = "/Volumes/LaCie/VenomCap/data-analysis/read-stats"
  # quiet = TRUE
  # overwrite = FALSE

  if (is.null(read.directory) == T){ stop("A directory of genome(s) is needed.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  # Gathers sample name data
  reads = list.files(read.directory, recursive = T, full.names = T)
  sample.names = list.dirs(read.directory, recursive = F, full.names = F)
  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  #Stats table prepare
  #Sets up data to collect
  header.data = c("sample", "total_lanes", "total_nucleotides", "total_gigabases", "total_read_pairs", "total_reads",
                  "read1_nucleotides", "read2_nucleotides", "read1_reads", "read2_reads")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, sample:=as.character(sample)]

  #Loops through each locus and does operations on them
  for (i in 1:length(sample.names)){

    #Reads in the contig file to check
    dir.create(paste0(output.directory, "/", sample.names[i]))
    #     sample.contigs = paste0(read.directory, "/", sample.names[i], "/assembled-contigs/", sample.names[i], "_orthologs.fa")

    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]

    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }

    sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    #Sets up data to collect
    s.header.data = c("sample", "total_nucleotides", "total_gigabases", "total_read_pairs", "total_reads",
                    "read1_nucleotides", "read2_nucleotides", "read1_reads", "read2_reads")
    sample.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.reads), ncol = length(s.header.data)))
    data.table::setnames(sample.data, s.header.data)
    sample.data[, sample:=as.character(sample)]

    for (j in 1:length(sample.reads)){

      lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]

      #Checks the Sample column in case already renamed
      if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }
      #Returns an error if reads are not found
      if (length(lane.reads) == 0 ){
        stop(sample.reads[j], " does not have any reads present for files ")
      } #end if statement

      #Copies ortholog file to act as a reference
      #       system(paste0("cp ", sample.contigs, " ", output.directory, "/", sample.names[i], "/sample-reference.fa"))

      r.reads1 = ShortRead::countFastq(lane.reads[1])
      r.reads2 = ShortRead::countFastq(lane.reads[2])

      data.table::set(sample.data, i = as.integer(j), j = match("sample", s.header.data), value = sample.names[i])
      data.table::set(sample.data, i = as.integer(j), j = match("total_nucleotides", s.header.data), value = r.reads1$nucleotides + r.reads2$nucleotides)
      data.table::set(sample.data, i = as.integer(j), j = match("total_gigabases", s.header.data), value = (r.reads1$nucleotides + r.reads2$nucleotides) / 1000000000 )
      data.table::set(sample.data, i = as.integer(j), j = match("total_read_pairs", s.header.data), value = r.reads1$records )
      data.table::set(sample.data, i = as.integer(j), j = match("total_reads", s.header.data), value = r.reads1$records+r.reads2$records )
      data.table::set(sample.data, i = as.integer(j), j = match("read1_nucleotides", s.header.data), value = r.reads1$nucleotides )
      data.table::set(sample.data, i = as.integer(j), j = match("read2_nucleotides", s.header.data), value = r.reads2$nucleotides )
      data.table::set(sample.data, i = as.integer(j), j = match("read1_reads", s.header.data), value = r.reads1$records )
      data.table::set(sample.data, i = as.integer(j), j = match("read2_reads", s.header.data), value = r.reads2$records )

    }#end j loop

    write.table(sample.data, file = paste0(output.directory, "/", sample.names[i], "/sample_reads_raw-data.txt"), sep = "\t", row.names = F)


    data.table::set(collect.data, i = as.integer(i), j = match("sample", header.data), value = sample.names[i])
    data.table::set(collect.data, i = as.integer(i), j = match("total_lanes", header.data), value = nrow(sample.data))
    data.table::set(collect.data, i = as.integer(i), j = match("total_nucleotides", header.data), value = sum(sample.data$total_nucleotides) )
    data.table::set(collect.data, i = as.integer(i), j = match("total_gigabases", header.data), value = sum(sample.data$total_gigabases) )
    data.table::set(collect.data, i = as.integer(i), j = match("total_read_pairs", header.data), value = sum(sample.data$total_read_pairs) )
    data.table::set(collect.data, i = as.integer(i), j = match("total_reads", header.data), value = sum(sample.data$total_reads) )
    data.table::set(collect.data, i = as.integer(i), j = match("read1_nucleotides", header.data), value = sum(sample.data$read1_nucleotides) )
    data.table::set(collect.data, i = as.integer(i), j = match("read2_nucleotides", header.data), value = sum(sample.data$read2_nucleotides) )
    data.table::set(collect.data, i = as.integer(i), j = match("read1_reads", header.data), value = sum(sample.data$read1_reads) )
    data.table::set(collect.data, i = as.integer(i), j = match("read2_reads", header.data), value = sum(sample.data$read2_reads) )

    print(paste0(sample.names[i], " Completed!"))

  }# end i loop

  write.table(collect.data, file = paste0(output.directory, "/sample-read-summary.txt"), sep = "\t", row.names = F)

  #### Print some textual summary here

}#end function

#END SCRIPT
