#' @title genomeCoverage
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

genomeCoverage = function(genome.directory = NULL,
                          read.directory = NULL,
                          output.directory = NULL,
                          threads = 1,
                          memory = 1,
                          overwrite = FALSE,
                          resume = TRUE,
                          quiet = TRUE,
                          samtools.path = NULL,
                          bwa.path = NULL,
                          picard.path = NULL) {


  #Read in basic genome info
  setwd("/Volumes/Rodents/Shrew_Genome")
  genome.file= "/Volumes/Rodents/Shrew_Genome/Draft_assemblies/HH664_spades.fa"
  read.directory = "/Volumes/Rodents/Shrew_Genome/processed-reads/Test"
  output.directory = "/Volumes/Rodents/Shrew_Genome/genomeStats/CovSim"
  threads = 4
  memory = 4
  samtools.path = "/Users/chutter/miniconda3/bin"
  bwa.path = "/usr/local/bin"
  picard.path = "/Users/chutter/miniconda3/bin"
  bbmap.path = "/usr/local/bin"
  resume = F
  overwrite = T
  max.depth = 50
  overwrite.reference = T
  number.gb = c(1, 2, 3, 4,5, 6, 7, 8, 9, 10)
  read.length = 150
  random.contig.count = 1000

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
  if (is.null(picard.path) == FALSE){
    b.string = unlist(strsplit(picard.path, ""))
    if (b.string[length(b.string)] != "/") {
      picard.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { picard.path = "" }

  #Same adds to bbmap path
  if (is.null(bbmap.path) == FALSE){
    b.string = unlist(strsplit(bbmap.path, ""))
    if (b.string[length(b.string)] != "/") {
      bbmap.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bbmap.path = "" }


  if (is.null(genome.file) == T){ stop("A directory of genome(s) is needed.") }

  #So I don't accidentally delete everything while testing resume
  if (resume == TRUE & overwrite == TRUE){
    overwrite = FALSE
    stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
  }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #genome.files = list.files(genome.directory)

  #Stats table prepare
  header.data = c("assembly", "n_contigs_total", "n_contigs_cover", "n_bp_total", "n_bp_cover",
                  "mean_depth", "median_depth",
                  "no_coverage", "above_1x", "above_2x", "above_3x", "above_4x", "above_5x", "above_10x", "above_20x")
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(number.gb), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, assembly:=as.character(assembly)]

  ### Make a reference function
  #Copies ortholog file to act as a reference
  ref.path = paste0(output.directory, "/ref_index")
  dir.create(ref.path)
  sample.name = gsub(".*/", "", genome.file)
  sample.name = gsub(".fa$", "", sample.name)
  system(paste0("cp ", genome.file, " ", ref.path, "/reference.fa"))

  #BWA indexes reads
  system(paste0(bwa.path, "bwa index -p ", ref.path, "/reference ", ref.path, "/reference.fa"),
         ignore.stderr = quiet, ignore.stdout = quiet)
  system(paste0(samtools.path, "samtools faidx ", ref.path, "/reference.fa"),
         ignore.stderr = quiet, ignore.stdout = quiet)
  system(paste0(picard.path, "picard -Xmx", memory, "G",
                " CreateSequenceDictionary REFERENCE=", ref.path, "/reference.fa",
                " OUTPUT=", ref.path, "/reference.dict USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

  #BWA mapping
  read.files = list.files(read.directory)
  read.1 = paste0(read.directory, "/", read.files[grep("READ1|_R1|_1_", read.files)])
  read.2 = paste0(read.directory, "/", read.files[grep("READ2|_R2|_2_", read.files)])

  #### Count total reads in file
  read.count = as.numeric(system(paste0("zcat < ", read.1, " | echo $((`wc -l`/4))"), intern = T))

  #Simulates coverage
  for (i in 1:length(number.gb)){

    #Gets the amount of reads to randomly replace
    number.reads = number.gb[i]*1000000000/(read.length * 2)
    sample.value = number.reads/read.count

    #Randomly draw out that many reads from the originals and save
    subsample.name = paste0(sample.name, "_", number.gb[i])
    suboutput = paste0(output.directory, "/", subsample.name)
    dir.create(paste0(output.directory, "/", subsample.name))

    #Uses bbmap reformat for this
    system(paste0(bbmap.path, "reformat.sh in1=", read.1, " in2=", read.2,
                  " out1=", suboutput, "/", subsample.name, "_READ1.fastq.gz",
                  " out2=", suboutput, "/", subsample.name, "_READ2.fastq.gz",
                  " samplerate=", sample.value, " overwrite=true"))

    #Subreads
    subread.1 = paste0(suboutput, "/", subsample.name, "_READ1.fastq.gz")
    subread.2 = paste0(suboutput, "/", subsample.name, "_READ2.fastq.gz")

    #Runs BWA mapping new read subset
    system(paste0(bwa.path, "bwa mem -M -E -0 -k 100 -w 4 -L 100",
                  " -t ", threads, " ", ref.path, "/reference ",
                  subread.1, " ", subread.2, " | ", samtools.path ,"samtools sort -@ ", threads,
                  " -o ", output.directory, "/mapped-reads.bam"),
           ignore.stderr = quiet, ignore.stdout = quiet)

    #Take mapped reads file and get coverage somehow from bbap?
    #system(paste0(bbmap.path, "pileup.sh in=", output.directory, "/mapped-reads.bam",
    #              " out=", output.directory, "/coverage-stats.txt overwrite=true"))

    # system(paste0(samtools.path, "samtools depth ", output.directory, "/mapped-reads.bam",
    #  " > ", output.directory, "/coverage-stats.txt"))

    system(paste0(samtools.path, "samtools index ", output.directory, "/mapped-reads.bam"),
           ignore.stderr = quiet, ignore.stdout = quiet)

    ### Make a temp bed file
    genome.data = Biostrings::readDNAStringSet(genome.file)
    genome.data = genome.data[Biostrings::width(genome.data) >= 500]
    rand.sampl = genome.data[sample(x = length(genome.data), size = random.contig.count, replace = FALSE)]
    rand.head = c("chrom", "chromStart", "chromEnd")
    save.bed = data.table::data.table(matrix(as.numeric(0), nrow = length(rand.sampl), ncol = length(rand.head)))
    data.table::setnames(save.bed, rand.head)
    save.bed[, chrom:=as.character(names(rand.sampl))]
    save.bed[, chromStart:=1]
    save.bed[, chromEnd:=Biostrings::width(rand.sampl)]
    write.table(save.bed, file = paste0(output.directory, "/random-target.bed"),
                row.names = F, quote = F, sep = "\t", col.names = F)

    system(paste0(samtools.path, "mosdepth --no-per-base --fast-mode --use-median ",
                  "--by ", output.directory, "/random-target.bed --threads ", threads,
                  " ", output.directory, "/", subsample.name, " ", output.directory, "/mapped-reads.bam"))

    #Reads in depth information
    depth.summary = data.table::fread(paste0(output.directory, "/", subsample.name, ".mosdepth.summary.txt"), header = T)
    depth.summary = depth.summary[grep("_region", depth.summary$chrom, invert = T),]

    #Saves the data
    data.table::set(collect.data, i = as.integer(i), j = match("assembly", header.data), value = subsample.name)
    data.table::set(collect.data, i = as.integer(i), j = match("n_contigs_total", header.data), value = length(rand.sampl) )
    data.table::set(collect.data, i = as.integer(i), j = match("n_contigs_cover", header.data), value = nrow(depth.summary) )
    data.table::set(collect.data, i = as.integer(i), j = match("n_bp_total", header.data), value = sum(Biostrings::width(rand.sampl)) )
    data.table::set(collect.data, i = as.integer(i), j = match("n_bp_cover", header.data), value = sum(depth.summary$bases) )
    data.table::set(collect.data, i = as.integer(i), j = match("mean_depth", header.data), value = mean(depth.summary$mean) )
    data.table::set(collect.data, i = as.integer(i), j = match("median_depth", header.data), value = median(depth.summary$mean) )

    #Depth plot at different thresholds
    header.depth = c("contig", "coverage", "proportion")
    region.summary = data.table::fread(paste0(output.directory, "/", subsample.name, ".mosdepth.region.dist.txt"), header = F)
    data.table::setnames(region.summary, header.depth)

    region.summary = region.summary[region.summary$coverage <= max.depth,]
    total.depth = aggregate(x = region.summary$proportion, by = list(region.summary$coverage), FUN = median)

    data.table::set(collect.data, i = as.integer(i), j = match("no_coverage", header.data), value = (total.depth$x[1]-total.depth$x[2])*100 )
    data.table::set(collect.data, i = as.integer(i), j = match("above_1x", header.data), value = total.depth$x[2]*100 )
    data.table::set(collect.data, i = as.integer(i), j = match("above_2x", header.data), value = total.depth$x[3]*100 )
    data.table::set(collect.data, i = as.integer(i), j = match("above_3x", header.data), value = total.depth$x[4]*100 )
    data.table::set(collect.data, i = as.integer(i), j = match("above_4x", header.data), value = total.depth$x[5]*100 )
    data.table::set(collect.data, i = as.integer(i), j = match("above_5x", header.data), value = total.depth$x[6]*100 )
    data.table::set(collect.data, i = as.integer(i), j = match("above_10x", header.data), value = total.depth$x[11]*100 )
    data.table::set(collect.data, i = as.integer(i), j = match("above_20x", header.data), value = total.depth$x[21]*100 )

  }#end i loop

  write.csv(collect.data, file = paste0(output.directory, "/", sample.name, "_coverage-summary.csv"), row.names = F)

  ##### Make plots
  mosdef.files = list.files(output.directory)
  mosdef.files = mosdef.files[grep("mosdepth.region", mosdef.files)]
  header.depth = c("contig", "coverage", "proportion")

  plot.data = c()
  for (i in 1:length(mosdef.files)){

    #Depth plot at different thresholds
    region.summary = data.table::fread(paste0(output.directory, "/", mosdef.files[i]), header = F)
    data.table::setnames(region.summary, header.depth)

    region.summary = region.summary[region.summary$coverage <= max.depth,]
    total.depth = aggregate(x = region.summary$proportion, by = list(region.summary$coverage), FUN = median)

    temp.data = data.frame(Analysis = gsub(".mosdepth.*", "", mosdef.files[i]),
                           Coverage = total.depth$Group.1,
                           Proportion = total.depth$x)

    plot.data = rbind(plot.data, temp.data)

  }#end i loop

  library(ggplot2)

  plot.data = plot.data[plot.data$Coverage <= 20,]

  ggplot(data = plot.data, aes(x=Coverage, y=Proportion, group=Analysis, color=Analysis)) +
    geom_line()


} #End function

