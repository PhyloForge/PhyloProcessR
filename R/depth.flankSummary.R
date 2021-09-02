#' @title flankDepthSummary
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

flankDepthSummary = function(depth.directory = NULL,
                             output.directory = "flank-depth-summary",
                             target.file = NULL,
                             sample.groups = NULL,
                             threads = 1,
                             memory = 1,
                             overwrite = FALSE) {

  #Debug
  setwd("/Volumes/Armored/FrogCap_Anura_Seqcap")
  depth.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Analyses/read-depth"
  sample.groups = "sample_probeset_group.csv"
  output.directory = "Analyses/flank-depth-summary"
  target.file = "Somehow all 3"
  overwrite = TRUE
  max.value = 40000000

  ### OLD Stuff
  mean.data = fread("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files/Ranoidea_marker_mean-depth.txt")
  median.data = fread("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files/Ranoidea_marker_median-depth.txt")
  align.dir = "/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Post_Processing_New/Final_Alignments"
  work.dir = paste0(align.dir,  "/Alignments_Ranoidea/all-markers_untrimmed")
  setwd(work.dir)
  locus.names = list.files(work.dir)

  library(ggplot2)
  parameters = theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
                                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(vjust = 1, hjust = 0.5, size = 28)) +
    theme(axis.title.x = element_text(size = 20, vjust=-0.5),
          axis.title.y = element_text(size = 20, vjust=2.2)) +
    theme(text = element_text(size=20),
          axis.text.x = element_text(vjust=1))


  #Quick checks
  if (is.null(depth.directory) == TRUE){ stop("Please provide input directory.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (file.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Read in sample data **** sample is run twice?!
  stat.files = list.files(depth.directory, recursive = T, full.names = T)
  if (is.null(sub.directory) == FALSE) {
    stat.files = stat.files[grep(paste0(sub.directory, "/"), stat.files)]
    sample.names = gsub(paste0("/", sub.directory, "/.*"), "", stat.files)
    sample.names = unique(gsub(paste0(depth.directory, "/"), "", sample.names))
  } else {
    sample.names = list.files(depth.directory, recursive = F, full.names = F)
  }

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  group.data = read.csv(sample.groups, header = T)
  group.data = group.data[group.data[,2] %in% sample.names,]
  group.names = unique(group.data[,1])

  ##################################### Sample Summary stats
  ###############################################################################################
  #Sets up data to collect
  header.spp = c("group", "sample", "mean_rpkm", "median_rpkm", "sd_rpkm", "min_rpkm", "max_rpkm")
  spp.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.spp)))
  data.table::setnames(spp.data, header.spp)
  spp.data[, sample:=as.character(sample)]
  spp.data[, group:=as.character(group)]

  for (j in 1:length(group.names)){

    #Runs through each sample
    #all.data = data.table::data.table()

    group.samples = group.data[group.data[,1] %in% group.names[j],][,2]
    bin.vals = vector("list", 100)
    for (i in 1:length(group.samples)) {
      #################################################
      ### Part A: prepare for loading and checks
      #################################################
      if (is.null(sub.directory) == TRUE){
        sample.data = as.matrix( data.table::fread(paste0(depth.directory, "/", group.samples[i], "/rpkm-binned-data.txt")) )
      } else {
        sample.data = as.matrix( data.table::fread(paste0(depth.directory, "/",group.samples[i], "/", sub.directory, "/rpkm-binned-data.txt")) )
      }

      sample.median = apply(sample.data, 1, median, na.rm=T)
      bin.vals = mapply(c, bin.vals, sample.median, SIMPLIFY=FALSE)

    } #end i loop

    #Gets barplot data.frame ready
    bin.data = data.frame(bin = seq(1:100), median = unlist(lapply(bin.vals, median)),
                         sd = unlist(lapply(bin.vals, sd)) )

    #ggplot rocks
    ggplot(bin.data, aes(x=bin, y=median)) + parameters +
      geom_bar(position=position_dodge(), stat="identity", size=.3, color = "#00FF7F", fill = "#00FF7F") +
      geom_errorbar(aes(ymin=median-sd, ymax=median+sd), size=.1, width=.1, position=position_dodge(.9), color = "black") +
      ggtitle(group.names[j]) + xlab("Contig (1 percent bins)") + ylab("Median Depth (RPKM)") + ylim(0, max.value)

    ggplot2::ggsave(paste0(output.directory, "/", group.names[j], "_bin-plot_rpkm-median.pdf"), width = 8, height = 6)

  }#end j loop

#To do start here?
  ### OLD Stuff
  mean.data = fread("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files/Ranoidea_marker_mean-depth.txt")
  median.data = fread("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files/Ranoidea_marker_median-depth.txt")
  bait.loci = scanFa(FaFile("/Users/chutter/Dropbox/Research/2_WIP/FrogCap_Pipeline/Source_Files/Final_Probesets/Old/Final_Ranoidea_Probe-Loci_Nov20.fa"))
  align.dir = "/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Post_Processing_New/Final_Alignments"
  work.dir = paste0(align.dir,  "/Alignments_Ranoidea/all-markers_untrimmed")
  setwd(work.dir)
  locus.names = list.files(work.dir)

  #Gets the subdata
  #setwd(sample.folders[6])
  #sample.names<-list.files(".", recursive = F)
  sample.sheet = read_excel("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_Mar132019.xlsx")
  sample.sheet = sample.sheet[sample.sheet$FrogCap_Paper == 1,]
  sample.names = sample.sheet[sample.sheet$Ranoidea == 1,]$Sample

  #Sets up database
  header.data<-c("Sample", locus.names)
  collect.data<-data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
  setnames(collect.data, header.data)
  collect.data[, Sample:=as.character(sample.names)]

  mean.exon.data = c()
  mean.intron.data = c()
  median.exon.data = c()
  median.intron.data = c()
  for (i in 1:length(locus.names)){

    #Read in alignment
    temp.name = gsub(".phy$", "", locus.names[i])
    align = readAAMultipleAlignment(file = paste0(work.dir,"/", locus.names[i]), format = "phylip")

    #Removes samples that too short individually
    write.temp<-strsplit(as.character(align), "")
    aligned.set<-as.matrix(as.DNAbin(write.temp) )
    len.temp<-as.character(as.list(aligned.set))
    len.loci<-lapply(len.temp, function (x) x[x != "-"])
    aln.len<-unlist(lapply(len.loci, function (x) length(x)))

    #Add in probe file
    ref.locus = bait.loci[names(bait.loci) == temp.name]

    if (length(ref.locus) == 0){
      print("Not found")
      next }

    names(ref.locus) = "Reference"
    alignment = run.mafft(align, add.contigs = ref.locus, algorithm = "add")

    #Finds the locus data
    lo.no = match(gsub(".phy$", "", locus.names[i]), colnames(mean.data))
    temp.mean = mean.data[,..lo.no]
    temp.mean$Bin = seq(1:100)
    lo.no = match(gsub(".phy$", "", locus.names[i]), colnames(median.data))
    temp.median = median.data[,..lo.no]
    temp.median$Bin = seq(1:100)

    #Trim alignment to find quantity of intron data
    #Removes the edge gaps
    ref.aligned<-as.character(alignment['Reference'])
    not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
    new.align<-strsplit(as.character(alignment), "")
    x<-as.matrix(as.DNAbin(new.align))
    trim.align = x[,-not.gaps]

    #Converts alignment back to DNASTringSet
    char.align<-as.list(data.frame(t(as.character(trim.align))))
    temp.align<-lapply(char.align, FUN = function(x) paste(x, collapse = ""))
    trimmed<-DNAStringSet(unlist(temp.align))

    #Removes samples that too short individually
    write.temp<-strsplit(as.character(trimmed), "")
    aligned.set<-as.matrix(as.DNAbin(write.temp) )
    len.temp<-as.character(as.list(aligned.set))
    len.loci<-lapply(len.temp, function (x) x[x != "-"])
    spp.len<-unlist(lapply(len.loci, function (x) length(x)))
    spp.len = spp.len[names(spp.len) != "Reference"]
    spp.len = spp.len / aln.len

    #Collects where teh intron starts
    set(collect.data, i = match(names(spp.len), sample.names),
        j = match(locus.names[i], header.data), value = as.numeric(spp.len))

    #Next, collect data to separately estimate intron and exon depth
    #locus name, exon start, exon end? then use with other stuff to get depth?
    #Or use same marker cutter as above but remove bases from the coverage sheet?
    #***** USE BIN DATA HERE

    start = round(mean(spp.len / 2) * 100)
    end = 100 - start

    #Collect mean data
    exon.bin.data = temp.mean[temp.mean$Bin %in% start:end,]
    intron.bin.data = temp.mean[!temp.mean$Bin %in% exon.bin.data$Bin,]
    exon.bin.data = exon.bin.data[,1]
    intron.bin.data = intron.bin.data[,1]

    #Collects where teh intron starts
    mean.exon.data = append(mean.exon.data, colMeans(exon.bin.data) )
    mean.intron.data = append(mean.intron.data, colMeans(intron.bin.data) )

    #collect median data
    exon.bin.data = temp.median[temp.median$Bin %in% start:end,]
    intron.bin.data = temp.median[!temp.median$Bin %in% exon.bin.data$Bin,]
    exon.bin.data = exon.bin.data[,1]
    intron.bin.data = intron.bin.data[,1]

    #Collects where teh intron starts
    median.exon.data = append(median.exon.data, colMeans(exon.bin.data))
    median.intron.data = append(median.intron.data, colMeans(intron.bin.data))
  } #end i loop

  #Saves the data
  save.temp.data = data.frame(Dataset = "Ranoidea", Type = "Exon",
                              Mean = mean(median.exon.data), Median = median(median.exon.data),
                              SD = sd(median.exon.data), min = min(median.exon.data), max = max(median.exon.data) )

  save.all.data = rbind(save.all.data, save.temp.data)

  save.temp.data = data.frame(Dataset = "Ranoidea", Type = "Intron",
                              Mean = mean(median.intron.data), Median = median(median.intron.data),
                              SD = sd(median.intron.data), min = min(median.intron.data), max = max(median.intron.data) )

  save.all.data = rbind(save.all.data, save.temp.data)

  setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
  write.table(collect.data, "Ranoidea_intron_boundaries.txt", row.names = F)

  write.ex = data.frame(Locus = names(median.exon.data), Exon = median.exon.data)
  write.in = data.frame(Locus = names(median.intron.data), Intron = median.intron.data)
  write.all = merge(write.ex, write.in, by = "Locus")
  write.table(write.all, "Ranoidea_depth-by-locus.txt", row.names = F)

  #Counts stuff
  save.data = c()
  for (i in 2:length(colnames(collect.data))){

    #temp data
    temp.data = data.frame(collect.data[,..i])

    #exclude zeros
    temp.data = temp.data[temp.data != 0,]
    if (length(temp.data) == 0){ next }

    save.temp = mean(temp.data)/2
    save.data = append(save.data, save.temp)
  } #end i

  #Mean intron percent each side
  mean(save.data)
  sd(save.data)
  min(save.data)
  max(save.data)





}#end function





