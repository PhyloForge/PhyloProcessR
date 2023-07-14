#' @title summary.sampleSpecificity
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param input.reads path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param decontamination.path directory of genomes contaminants to scan samples
#'
#' @param samtools.path system path to samtools in case it can't be found
#'
#' @param bwa.path system path to bwa in case it can't be found
#'
#' @param threads number of computation processing threads
#'
#' @param mem amount of system memory to use
#'
#' @param resume TRUE to skip samples already completed
#'
#' @param overwrite TRUE to overwrite a folder of samples with output.dir
#'
#' @param quiet TRUE to supress screen output
#'
#' @return a new directory of adaptor trimmed reads and a summary of the trimming in logs/
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

analysis.sampleSpecificity = function(read.directory = NULL,
                                      target.file = NULL,
                                      output.directory = "sample-specificity",
                                      threads = 1,
                                      memory = 1,
                                      overwrite = FALSE,
                                      quiet = FALSE,
                                      bwa.path = NULL,
                                      gatk4.path = NULL,
                                      samtools.path = NULL) {

  # read.directory = "/Users/chutter/Dropbox/VenomCap_test_data/cleaned-reads"
  # target.file = "/Users/chutter/Dropbox/VenomCap_test_data/venom_loci_updated_Mar12_cdhit95_duplicate_exons_renamed_Feb2023_FINAL.fa"
  # output.directory = "/Users/chutter/Dropbox/VenomCap_test_data/sample-specificity"
  
  # samtools.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin/"
  # bwa.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin/"
  # gatk4.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin/"
  
  # quiet = TRUE
  # overwrite = FALSE
  # threads = 6
  # memory = 6

  ##### Program path check
  ####################################################################
  #Same adds to bbmap path
  if (is.null(samtools.path) == FALSE){
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { samtools.path = "" }

  #Same adds to bbmap path
  if (is.null(bwa.path) == FALSE){
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bwa.path = "" }

  #Same adds to bbmap path
  if (is.null(gatk4.path) == FALSE){
    b.string = unlist(strsplit(gatk4.path, ""))
    if (b.string[length(b.string)] != "/") {
      gatk4.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { gatk4.path = "" }
  ####################################################################

  #Quick checks
  if (is.null(read.directory) == TRUE){ stop("Please provide input directory.") }

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

  #reads in target markers
  target.markers = Biostrings::readDNAStringSet(target.file)

  #Specificity refers to the percentage of cleaned reads that can be mapped back to the target markers.

  #Sets up data to collect
  header.all = c("sample", "mapped_reads", "unmapped_mates",  "total_unmapped_reads",
                  "total_read_pairs",  "median_marker_rpkm", "mean_marker_rpkm", "sample_specificity")
  all.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.all)))
  data.table::setnames(all.data, header.all)
  all.data[, sample:=as.character(sample)]

  #Loops through each locus and does operations on them
  for (i in 1:length(sample.names)){

    #Reads in the contig file to check
    dir.create(paste0(output.directory, "/", sample.names[i]))
    #     sample.contigs = paste0(read.directory, "/", sample.names[i], "/assembled-contigs/", sample.names[i], "_orthologs.fa")

    #Sets up data to collect
    header.data = c("sample", "marker", "mapped_reads", "unmapped_mates",
                    "marker_total_reads", "total_unmapped", "marker_rpkm")
    collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(target.markers), ncol = length(header.data)))
    data.table::setnames(collect.data, header.data)
    collect.data[, sample:=as.character(sample)]
    collect.data[, marker:=as.character(marker)]

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

      system(paste0("cp ", target.file, " ", output.directory, "/", sample.names[i], "/sample-reference.fa"))
      system(paste0(bwa.path, "bwa index ", output.directory, "/", sample.names[i], "/sample-reference.fa"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      #Creates a bam alignment file of reads mapped to reference
      system(paste0(bwa.path, "bwa mem -M -t ", threads, " ",
                    output.directory, "/", sample.names[i], "/sample-reference.fa ",
                    lane.reads[1], " ", lane.reads[2],
                    " | ", samtools.path, "samtools sort -@", threads, " -O BAM",
                    " -o ", output.directory, "/", sample.names[i], "/paired.bam  -"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
                    " MarkDuplicates -I ", output.directory, "/", sample.names[i], "/paired.bam",
                    " -O ", output.directory, "/", sample.names[i], "/dedup_paired.bam",
                    " -M ", output.directory, "/", sample.names[i], "/metrics.txt",
                    " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      system(paste0("mv ", output.directory, "/", sample.names[i], "/metrics.txt ",
                    output.directory, "/", sample.names[i], "/duplication_stats.txt"))

      #HERe again
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
                    " BuildBamIndex -I ", output.directory, "/", sample.names[i], "/paired.bam",
                    " -O ", output.directory, "/", sample.names[i], "/paired.bai",
                    " -USE_JDK_DEFLATER true -USE_JDK_INFLATER true"),
             ignore.stdout = quiet, ignore.stderr = quiet)

      #Gets the genome coverage using bedtools
      system(paste0(samtools.path, "samtools idxstats -@", threads, " ",
                    output.directory, "/", sample.names[i], "/paired.bam > ",
                    output.directory, "/", sample.names[i], "/samtools_idxstats.txt"))

      #system(paste0("rm ", output.directory, "/", sample.names[i], "/paired.bam"))

      #IDX data
      i.headers = c("marker", "length", "mapped_reads", "unmapped_mates")
      idx.stats = data.table::fread(file = paste0(output.directory, "/", sample.names[i], "/samtools_idxstats.txt"))
      data.table::setnames(idx.stats, i.headers)
      idx.stats = idx.stats[idx.stats$marker != "*",]

      #Gets the number of reads that mapped to the reference, also total reads
      #Total number of un-mapped reads
      mapped.all = sum(idx.stats$mapped_reads)
      mapped.all = as.numeric(system(paste0(samtools.path, "samtools view -c -F 4 ",
                                             output.directory, "/", sample.names[i], "/paired.bam"), intern = T))
      unmapped.all = as.numeric(system(paste0(samtools.path, "samtools view -c -f 4 ",
                                              output.directory, "/", sample.names[i], "/paired.bam"), intern = T))

      #Goes through each locus and calculates stats
      locus.names = unique(idx.stats$marker)
      locus.rpkm = idx.stats$mapped_reads / (idx.stats$length/1000 * mapped.all/1000000)

      #Starter stats
      data.table::set(collect.data, j = match("sample", header.data), value = sample.names[i] )
      data.table::set(collect.data, j = match("marker", header.data), value = locus.names )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("mapped_reads", header.data), value = idx.stats$mapped_reads )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("unmapped_mates", header.data), value = idx.stats$unmapped_mates )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("marker_total_reads", header.data), value = (idx.stats$unmapped_mates+idx.stats$mapped_reads) )
      data.table::set(collect.data, j = match("total_unmapped", header.data), value = unmapped.all )
      data.table::set(collect.data, i = match(locus.names, idx.stats$marker), j = match("marker_rpkm", header.data), value = locus.rpkm )

    } #end j loop

    write.table(collect.data, file = paste0(output.directory, "/", sample.names[i], "/sample_specificity_raw-data.txt"), sep = "\t", row.names = F)

    #Collects the overall summary data
    collect.data = collect.data[collect.data$mapped_reads != 0,]
    data.table::set(all.data, i = as.integer(i), j = match("sample", header.all), value = sample.names[i] )
    data.table::set(all.data, i = as.integer(i), j = match("mapped_reads", header.all), value = sum(collect.data$mapped_reads) )
    data.table::set(all.data, i = as.integer(i), j = match("unmapped_mates", header.all), value = sum(collect.data$unmapped_mates) )
    data.table::set(all.data, i = as.integer(i), j = match("total_unmapped_reads", header.all), value = unmapped.all )
    data.table::set(all.data, i = as.integer(i), j = match("total_read_pairs", header.all), value = sum(collect.data$marker_total_reads)+unmapped.all )
    data.table::set(all.data, i = as.integer(i), j = match("mean_marker_rpkm", header.all), value = mean(collect.data$marker_rpkm) )
    data.table::set(all.data, i = as.integer(i), j = match("median_marker_rpkm", header.all), value = median(collect.data$marker_rpkm) )
    data.table::set(all.data, i = as.integer(i), j = match("sample_specificity", header.all), value = sum(collect.data$marker_total_reads)/(sum(collect.data$marker_total_reads)+unmapped.all) )

    print(paste0(sample.names[i], " Completed!"))

  }# end i loop

  write.table(all.data, file = paste0(output.directory, "/sample-specificity_summary.txt"), sep = "\t", row.names = F)

  #### Print some textual summary here

}# end function

### END SCRIPT
