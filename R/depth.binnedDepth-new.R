#' @title summarizeDepth
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

# binnedDepth = function(read.directory = NULL,
#                        sub.directory = NULL,
#                        output.name = "binned-depth",
#                        number.bins = 100,
#                        sample.groups = NULL,
#                        threads = 1,
#                        mem = 1,
#                        overwrite = FALSE) {
#
#   #Debug
#   setwd("/Volumes/Armored/FrogCap_Anura_Seqcap/Analyses")
#   read.directory = "/Volumes/Armored/FrogCap_Anura_Seqcap/Processed_Samples"
#   sub.directory = "cleaned-reads-snp"
#   output.name = "binned-depth"
#   overwrite = FALSE
#   mosdepth.path = "/Users/chutter/miniconda3/bin"
#   samtools.path = "/Users/chutter/miniconda3/bin"
#   bwa.path = "/usr/local/bin"
#   picard.path = "/Users/chutter/miniconda3/bin"
#   threads = 4
#   memory = 4
#   quiet = TRUE
#   number.bins = 100
#
#
#   ##### Program path check
#   ####################################################################
#   #Same adds to bbmap path
#   if (is.null(samtools.path) == FALSE){
#     b.string = unlist(strsplit(samtools.path, ""))
#     if (b.string[length(b.string)] != "/") {
#       samtools.path = paste0(append(b.string, "/"), collapse = "")
#     }#end if
#   } else { samtools.path = "" }
#
#   #Same adds to bbmap path
#   if (is.null(bwa.path) == FALSE){
#     b.string = unlist(strsplit(bwa.path, ""))
#     if (b.string[length(b.string)] != "/") {
#       bwa.path = paste0(append(b.string, "/"), collapse = "")
#     }#end if
#   } else { bwa.path = "" }
#
#   #Same adds to bbmap path
#   if (is.null(picard.path) == FALSE){
#     b.string = unlist(strsplit(picard.path, ""))
#     if (b.string[length(b.string)] != "/") {
#       picard.path = paste0(append(b.string, "/"), collapse = "")
#     }#end if
#   } else { picard.path = "" }
#   ####################################################################
#
#   #Quick checks
#   if (is.null(read.directory) == TRUE){ stop("Please provide input directory.") }
#
#   #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
#   if (file.exists(output.name) == F){ dir.create(output.name) } else {
#     if (overwrite == TRUE){
#       system(paste0("rm ", output.name))
#       dir.create(output.name)
#     }
#   }#end else
#
#   #Read in sample data **** sample is run twice?!
#   reads = list.files(read.directory, recursive = T, full.names = T)
#   if (is.null(sub.directory) != TRUE) {
#     reads = reads[grep(paste0(sub.directory, "/"), reads)]
#     sample.names = gsub(paste0("/", sub.directory, "/.*"), "", reads)
#     sample.names = unique(gsub(paste0(read.directory, "/"), "", sample.names))
#   } else {
#     sample.names = list.files(read.directory, recursive = F, full.names = F)
#   }
#
#   if (length(sample.names) == 0){ return("no samples remain to analyze.") }
#
#   #Sets up data to collect
#   header.data = c("Sample","numberLoci", "meanDepth", "medianDepth", "sdDepth", "minDepth", "maxDepth")
#   collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
#   data.table::setnames(collect.data, header.data)
#   collect.data[, Sample:=as.character(Sample)]
#
#   #Sets up bin data
#   bin.collect = data.table::data.table(matrix(as.numeric(0), nrow = 100, ncol = length(sample.names)))
#   data.table::setnames(bin.collect, gsub(".*\\/", "", sample.names))
#
#   #Loops through each locus and does operations on them
#   for (i in 1:length(sample.names)){
#
#     #Reads in the contig file to check
#     dir.create(paste0(output.name, "/", sample.names[i]))
#     sample.contigs = paste0(read.directory, "/", sample.names[i], "/assembled-contigs/", sample.names[i], "_orthologs.fa")
#
#     #################################################
#     ### Part A: prepare for loading and checks
#     #################################################
#     sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]
#
#     #Checks the Sample column in case already renamed
#     if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }
#
#     sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))
#
#     #Returns an error if reads are not found
#     if (length(sample.reads) == 0 ){
#       stop(sample.names[i], " does not have any reads present for files ")
#     } #end if statement
#
#     for (j in 1:length(sample.reads)){
#
#       lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]
#
#       #Checks the Sample column in case already renamed
#       if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }
#       #Returns an error if reads are not found
#       if (length(lane.reads) == 0 ){
#         stop(sample.reads[j], " does not have any reads present for files ")
#       } #end if statement
#
#       #Copies ortholog file to act as a reference
#       system(paste0("cp ", sample.contigs, " ", output.name, "/", sample.names[i], "/sample-reference.fa"))
#       system(paste0(bwa.path, "bwa index ", output.name, "/", sample.names[i], "/sample-reference.fa"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       #Creates a bam alignment file of reads mapped to reference
#       system(paste0(bwa.path, "bwa mem -M -E -0 -k 100 -w 4 -L 100 -t ", threads, " ",
#                     output.name, "/", sample.names[i], "/sample-reference.fa ",
#                     lane.reads[1], " ", lane.reads[2],
#                     " | ", samtools.path, "samtools sort -@", threads, " -O BAM",
#                     " -o ", output.name, "/", sample.names[i], "/paired.bam  -"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       system(paste0(picard.path, "picard -Xmx", memory, "g",
#                     " MarkDuplicates INPUT=", output.name, "/", sample.names[i], "/paired.bam",
#                     " OUTPUT=", output.name, "/", sample.names[i], "/dedup_paired.bam",
#                     " METRICS_FILE=", output.name, "/", sample.names[i], "/metrics.txt",
#                     " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       system(paste0("mv ", output.name, "/", sample.names[i], "/metrics.txt ",
#                     output.name, "/", sample.names[i], "/duplication_stats.txt"))
#
#       #HERe again
#       system(paste0(picard.path, "picard -Xmx", memory, "g",
#                     " BuildBamIndex INPUT=", output.name, "/", sample.names[i], "/dedup_paired.bam",
#                     " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"),
#              ignore.stdout = quiet, ignore.stderr = quiet)
#
#       # ALternative if it doesn't work
#       #system(paste0(samtools.path, "samtools index ", output.directory, "/mapped-reads.bam"),
#       #       ignore.stderr = quiet, ignore.stdout = quiet)
#
#       ### Create bed file
#       contig.data = Biostrings::readDNAStringSet(filepath = sample.contigs, format = "fasta")
#       contig.ranges = data.frame(locus = names(contig.data),
#                                  size = Biostrings::width(contig.data),
#                                  bin_width = floor(Biostrings::width(contig.data)/number.bins))
#
#       contig.ranges = contig.ranges[contig.ranges$size >= 200,]
#
#       seqFun.start = function(a) {seq(from = 1, by = a, length.out = number.bins) }
#       seqFun.end = function(a) {seq(from = a, by = a, length.out = number.bins) }
#
#       bin.list.s = apply(contig.ranges[,2:3], 1, function(x) seqFun.start(x[2]))
#       names(bin.list.s) = contig.ranges$locus
#
#       bin.list.e = apply(contig.ranges[,2:3], 1, function(x) seqFun.end(x[2]))
#       names(bin.list.e) = contig.ranges$locus
#
#       # seqFun.start = function(a, b) {seq(from = 1, to = a-b, by = b) }
#       # seqFun.end = function(a, b) {seq(from = b, to = a, by = b) }
#       #
#       # bin.list.s = apply(contig.ranges[,2:3], 1, function(x) seqFun.start(x[1], x[2]))
#       # names(bin.list.s) = contig.ranges$locus
#       #
#       # bin.list.e = apply(contig.ranges[,2:3], 1, function(x) seqFun.end(x[1], x[2]))
#       # names(bin.list.e) = contig.ranges$locus
#       #
#       # temp.s = unlist(bin.list.s)
#
#       bed.file = data.frame(chrom = names(bin.list.e),
#                             chromStart = unlist(bin.list.s),
#                             chromEnd = unlist(bin.list.e) )
#
#       header.data = c("chrom", "chromStart", "chromEnd")
#       bed.data = data.table::data.table(matrix(as.numeric(0), nrow = length(unlist(bin.list.s)), ncol = length(header.data)))
#       data.table::setnames(bed.data, header.data)
#       bed.data[, chromStart := as.numeric(unlist(bin.list.s))]
#       bed.data[, chromEnd := as.numeric(unlist(bin.list.e))]
#       bed.data[, chrom := as.character(names(unlist(bin.list.e, use.names = T)))]
#       bed.data = bed.data[is.na(bed.data$chrom) != TRUE,]
#       ### Save bed file
#       write.table(x = bed.data)
#
#       #Change last one to length
#
#
#
#       bed.file = data.frame(chrom = names(bin.list), )
#
#       for (k in 1:nrow(contig.ranges)){
#
#
#
#
#         bin.ranges = seq(1:contig.ranges$size, by = contig.ranges$bin_width)
#
#
#         group.tags = cut(c(1, contig.ranges$size[k]),
#                          breaks = 100, include.lowest = T)
#
#         cut(contig.ranges$size, 100)
#
#         binr::bins(seq(1:contig.ranges$size[k]),
#                    target.bins = 100,
#                    max.breaks = NA,
#                    exact.groups = T,
#                    verbose = F,
#                    errthresh = 0.1,
#                    minpts = NA)
#
#         #1% bin fill in with mean
#         bin.size = floor(nrow(locus.stats)/number.bins)
#         start = 1
#         stop = bin.size
#         for (k in 1:number.bins){
#           loc.bin = locus.stats[start:stop,]
#           data.table::set(bin.data, i = as.integer(k), j = as.integer(j), value = median(loc.bin$mapped_reads) )
#           #new counters
#           start = start + bin.size
#           stop = stop + bin.size
#         }#end k loop
#
#       }
#
#
#
#
#       group.tags = cut(contig.ranges$size,
#                        breaks = contig.ranges$bin_width, include.lowest = TRUE)
#
#
#
#       group_tags <- cut(v$MeanEducation,
#                         breaks=breaks,
#                         include.lowest=TRUE,
#                         right=FALSE,
#                         labels=tags)
#
#
#       system(paste0(mosdepth.path, "mosdepth -n --fast-mode --use-median --by ", window.size,
#                     " --threads ", threads, " ", sample.names[i], " ",
#                     output.name, "/", sample.names[i], "/mapped-reads.bam"),
#              ignore.stderr = quiet, ignore.stdout = quiet)
#
#       #Reads in depth information
#       depth.summary = data.table::fread(paste0( output.name, "/", sample.names[i], "/",
#                                                 sample.names[i], ".mosdepth.summary.txt"), header = T)
#       depth.summary = depth.summary[grep("_region", depth.summary$chrom, invert = T),]
#
#       #Saves the data
#       i = 1
#       data.table::set(collect.data, i = as.integer(i), j = match("assembly", header.data), value = sample.name)
#       data.table::set(collect.data, i = as.integer(i), j = match("n_contigs", header.data), value = nrow(depth.summary) )
#       data.table::set(collect.data, i = as.integer(i), j = match("total_bp", header.data), value = sum(depth.summary$bases) )
#       data.table::set(collect.data, i = as.integer(i), j = match("mean_depth", header.data), value = mean(depth.summary$mean) )
#       data.table::set(collect.data, i = as.integer(i), j = match("median_depth", header.data), value = median(depth.summary$mean) )
#
#       #Depth plot at different thresholds
#       header.depth = c("contig", "coverage", "proportion")
#       depth.summary = data.table::fread(paste0(sample.name, ".mosdepth.global.dist.txt"), header = F)
#       data.table::setnames(depth.summary, header.depth)
#
#       depth.summary = depth.summary[depth.summary$coverage <= max.depth,]
#       total.depth = aggregate(x = depth.summary$proportion, by = list(depth.summary$coverage), FUN = median)
#
#       data.table::set(collect.data, i = as.integer(i), j = match("no_coverage", header.data), value = (total.depth$x[1]-total.depth$x[2])*100 )
#       data.table::set(collect.data, i = as.integer(i), j = match("above_1x", header.data), value = total.depth$x[2]*100 )
#       data.table::set(collect.data, i = as.integer(i), j = match("above_2x", header.data), value = total.depth$x[3]*100 )
#       data.table::set(collect.data, i = as.integer(i), j = match("above_3x", header.data), value = total.depth$x[4]*100 )
#       data.table::set(collect.data, i = as.integer(i), j = match("above_4x", header.data), value = total.depth$x[5]*100 )
#       data.table::set(collect.data, i = as.integer(i), j = match("above_5x", header.data), value = total.depth$x[6]*100 )
#       data.table::set(collect.data, i = as.integer(i), j = match("above_10x", header.data), value = total.depth$x[11]*100 )
#       data.table::set(collect.data, i = as.integer(i), j = match("above_20x", header.data), value = total.depth$x[21]*100 )
#
#       write.csv(collect.data, file = paste0(output.directory, "/", sample.name, "_coverage-summary.csv"), row.names = F)
#
#
#   }
#
#
#
#
#
#   group.data = read.csv(sample.groups, header = T)
#
#   base.headers = c("locus", "position", "mapped_reads")
#
#
#
#   for (j in 1:length(sample.names)){
#
#     #################################################
#     ### Part A: prepare for loading and checks
#     #################################################
#     #Load in perbase depth data
#     if (is.null(sub.directory) == TRUE){
#       sample.data = data.table::fread(paste0(depth.directory, "/", sample.names[i], "/samtools_perbase_depth.txt"))
#     } else {
#       sample.data = data.table::fread(paste0(depth.directory, "/",sample.names[i], "/", sub.directory, "/samtools_perbase_depth.txt"))
#     }
#     #give headers
#     data.table::setnames(sample.data, base.headers)
#     sample.data$locus = gsub("_\\|_.*", "", sample.data$locus)
#     locus.names = unique(sample.data$locus)
#
#     #Sets up database
#     header.data = c("bin", locus.names)
#     bin.data = data.table::data.table(matrix(as.numeric(0), nrow = number.bins, ncol = length(header.data)))
#     data.table::setnames(bin.data, header.data)
#     bin.data[, bin:=as.character(1:number.bins)]
#
#     #For each locus, caculate amount in bins
#     for (j in 1:length(locus.names)){
#
#       locus.stats = sample.data[sample.data$locus %in% locus.names[j]]
#
#       #1% bin fill in with mean
#       bin.size = floor(nrow(locus.stats)/number.bins)
#       start = 1
#       stop = bin.size
#       for (k in 1:number.bins){
#         loc.bin = locus.stats[start:stop,]
#         data.table::set(bin.data, i = as.integer(k), j = as.integer(j), value = median(loc.bin$mapped_reads) )
#         #new counters
#         start = start + bin.size
#         stop = stop + bin.size
#       }#end k loop
#
#   }#end j loop
#
#   #Save the two datasets inside the depth_analysis folder
#   write.table(bin.data, file = paste0(depth.directory, sample.names[i], "/binned_median_data.txt"),
#               sep = "\t", row.names = F)
#
#
#     #Gets group subset
#     group.stats = all.data[all.data$group %in% group.names[k],]
#     group.stats$locus = gsub("_\\|_.*", "", group.stats$locus)
#     locus.names = unique(group.stats$locus)
#
#     header.spp = c("group", "locus", "mean_length", "mean_rpkm", "median_rpkm", "sd_rpkm", "min_rpkm", "max_rpkm")
#     spp.data = data.table::data.table(matrix(as.numeric(0), nrow = length(locus.names), ncol = length(header.spp)))
#     data.table::setnames(spp.data, header.spp)
#     spp.data[, group:=as.character(group)]
#     spp.data[, locus:=as.character(locus)]
#
#
#     #Sets up database
#     header.data = c("Bin", locus.names)
#     locus.data = data.table(matrix(as.numeric(0), nrow = 100, ncol = length(header.data)))
#     setnames(locus.data, header.data)
#     locus.data[, Bin:=as.character(1:100)]
#
#     #Separate datasets
#     mean.data = locus.data
#     median.data = locus.data
#
#     #Sets up the list within a list
#     locus.list = vector("list", length(locus.names))
#     names(locus.list) = locus.names
#     the.bins = vector("list", 100)
#     for (i in 1:length(locus.list)){ locus.list[[i]] = the.bins }
#
#   #### end to do
#
#     locus.names = gsub(".phy$", "", list.files(work.dir))
#
#     #Gets the subdata
#     scale.data<-spp.data[spp.data$ProbeSet %in% "Hyloidea",]
#     sample.names<-scale.data$Sample
#
#     #Get ranoidea summary stats
#     hyl.data = spp.data[spp.data$Sample %in% sample.names,]
#
#     #Sets up database
#     header.data<-c("Bin", locus.names)
#     locus.data<-data.table(matrix(as.numeric(0), nrow = 100, ncol = length(header.data)))
#     setnames(locus.data, header.data)
#     locus.data[, Bin:=as.character(1:100)]
#
#     #Separate datasets
#     mean.data = locus.data
#     median.data = locus.data
#
#     #Sets up the list within a list
#     locus.list = vector("list", length(locus.names))
#     names(locus.list) = locus.names
#     the.bins = vector("list", 100)
#     for (i in 1:length(locus.list)){ locus.list[[i]] = the.bins }
#
#     #Loops through each samples binned median depth
#     for (j in 1:length(sample.names)){
#
#       #Go to file location
#       setwd(paste0(sample.folders[5], "/", sample.names[j], "/coverage-analysis"))
#       raw.cov<-as.matrix(fread("binned_median_data.txt"))
#
#       #Find a new way to match all the loci at once
#       for (k in 1:ncol(raw.cov)){
#         list.pos = match(colnames(raw.cov)[k], names(locus.list))
#         if (is.na(list.pos) == T){ next }
#         locus.list[[list.pos]] = mapply(c, locus.list[[list.pos]], raw.cov[,k], SIMPLIFY=FALSE)
#       }#end k loop
#
#     }#end j loop
#
#     #Summarizes data
#     for (i in 1:length(locus.list)){
#       #Gets barplot data.frame ready
#       bin.summ<-data.frame(Bin = seq(1:100), Mean = unlist(lapply(locus.list[[i]], mean)),
#                            Median = unlist(lapply(locus.list[[i]], median)), SD = unlist(lapply(locus.list[[i]], sd)) )
#
#       #Summarize the locus level
#       set(mean.data, i = 1:100, j = match(locus.names[i], header.data), value = bin.summ$Mean )
#       set(median.data, i = 1:100, j = match(locus.names[i], header.data), value = bin.summ$Median )
#     }#end i loop
#
#     #Saves the data
#     setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
#     write.table(mean.data, "Hyloidea_marker_mean-depth.txt", row.names = F)
#     write.table(median.data, "Hyloidea_marker_median-depth.txt", row.names = F)
#
#     #Makes the figure part
#     bin.vals<-vector("list", 100)
#     for (j in 1:length(sample.names)){
#       #Go to file location
#       setwd(paste0(sample.folders[5], "/", sample.names[j], "/coverage-analysis"))
#       raw.cov<-as.matrix(fread("binned_median_data.txt"))
#
#       bin.temp<-rowMedians(raw.cov, na.rm = T)
#       #bin.vals<-append(bin.vals, bin.temp)
#       bin.vals<-mapply(c, bin.vals, bin.temp, SIMPLIFY=FALSE)
#     }#end j loop
#
#     #Gets barplot data.frame ready
#     bin.data<-data.frame(Bin = seq(1:100), Median = unlist(lapply(bin.vals, median)),
#                          SD = unlist(lapply(bin.vals, sd)) )
#
#     #ggplot rocks
#     p1<-ggplot(bin.data, aes(x=Bin, y=Median)) + parameters +
#       geom_bar(position=position_dodge(), stat="identity", size=.3, color = "#00FF7F", fill = "#00FF7F") +
#       geom_errorbar(aes(ymin=Median-SD, ymax=Median+SD), size=.1, width=.1, position=position_dodge(.9), color = "black") +
#       ggtitle("Hyloidea") + xlab("Contig (1 percent bins)") + ylab("Median Depth of Coverage (X-bp)") + ylim(0,80)
#
#     p1
#     setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
#     ggsave("R_Plots/3_Hyloidea_Coverage.pdf", width = 6, height = 6)
#
#
#
#
#
#
#
#
#
#
# # RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
# #  numReads        - number of reads mapped to a gene sequence
# #  geneLength      - length of the gene sequence
# #  totalNumReads   - total number of mapped reads of a sample
#
#
#   #Creates the summary log
#   summary.data =  data.frame(Sample = as.character(),
#                              Read1_Count = as.numeric(),
#                              Read2_Count = as.numeric(),
#                              Total_Reads = as.numeric(),
#                              Read_Pairs = as.numeric(),
#                              Read_Length = as.numeric(),
#                              MegaBasePairs = as.numeric(),
#                              Reads_Per_Million = as.numeric())
#
#   #Runs through each sample
#   for (i in 1:length(sample.names)) {
#     #################################################
#     ### Part A: prepare for loading and checks
#     #################################################
#     if (is.null(sub.directory) == TRUE){
#       spp.data = data.table::fread(paste0(depth.directory, "/", sample.names[i], "/species_summary_data.txt"))
#     } else {
#       spp.data = data.table::fread(paste0(depth.directory, "/",sample.names[i], "/", sub.directory, "/species_summary_data.txt"))
#     }
#
#     locus.names = spp.data$locus
#
#     #Find a new way to match all the loci at once
#     for (k in 1:ncol(raw.cov)){
#       list.pos = match(colnames(raw.cov)[k], names(locus.list))
#       if (is.na(list.pos) == T){ next }
#       locus.list[[list.pos]] = mapply(c, locus.list[[list.pos]], raw.cov[,k], SIMPLIFY=FALSE)
#     }#end k loop
#
#   }#end j loop
#
#   #Summarizes data
#   for (i in 1:length(locus.list)){
#     #Gets barplot data.frame ready
#     bin.summ<-data.frame(Bin = seq(1:100), Mean = unlist(lapply(locus.list[[i]], mean)),
#                          Median = unlist(lapply(locus.list[[i]], median)), SD = unlist(lapply(locus.list[[i]], sd)) )
#
#     #Summarize the locus level
#     set(mean.data, i = 1:100, j = match(locus.names[i], header.data), value = bin.summ$Mean )
#     set(median.data, i = 1:100, j = match(locus.names[i], header.data), value = bin.summ$Median )
#   }#end i loop
#
#
#
#     sample.stats = stat.files[grep(pattern = paste0(sample.names[i], "_"), x = stat.files)]
#
#     #Checks the Sample column in case already renamed
#     if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = sample.names[i], x = reads)] }
#
#     sample.reads = unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))
#
#     #Returns an error if reads are not found
#     if (length(sample.reads) == 0 ){
#       stop(sample.names[i], " does not have any reads present for files ")
#     } #end if statement
#
#     for (j in 1:length(sample.reads)){
#
#       lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]
#
#       #Checks the Sample column in case already renamed
#       if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = sample.reads[j], x = reads)] }
#       #Returns an error if reads are not found
#       if (length(lane.reads) == 0 ){
#         stop(sample.reads[j], " does not have any reads present for files ")
#       } #end if statement
#
#       #Gathers stats on initial data
#       read1.count = as.numeric(system(paste0("zcat < ", lane.reads[1], " | echo $((`wc -l`/4))"), intern = T))
#       read2.count = as.numeric(system(paste0("zcat < ", lane.reads[2], " | echo $((`wc -l`/4))"), intern = T))
#       scale.factor = (read1.count + read2.count) / 1000000
#       if (read1.count == read2.count){ read.pairs = read1.count } else { read.pairs = "PairsUnequal"}
#
#       temp.remove = data.frame(Sample = sample.names[i],
#                                Read1_Count = read1.count,
#                                Read2_Count = read2.count,
#                                Total_Reads = read1.count + read2.count,
#                                Read_Pairs = read.pairs,
#                                Read_Length = read.length,
#                                MegaBasePairs = (read.length * read.pairs * 2)/1000000,
#                                Reads_Per_Million = scale.factor)
#
#       summary.data = rbind(summary.data, temp.remove)
#
#
#     }#end sample j loop
#
#     print(paste0(sample.names[i], " Completed fastq counting!"))
#
#   }#end sample i loop
#
#   write.csv(summary.data, file = paste0(output.name, ".csv"), row.names = FALSE)
# }
#
