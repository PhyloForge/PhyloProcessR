#' @title summary.geneticDistance
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

geneticDistance = function(alignment.directory = NULL,
                        target.file = NULL,
                        output.directory = "genetic-distance",
                        threads = 1,
                        memory = 1,
                        overwrite = FALSE,
                        quiet = FALSE,
                        mafft.path = NULL) {

      # alignment.directory = "/Users/chutter/Dropbox/VenomCap_test_data/untrimmed_all-markers"
      # target.file = "/Users/chutter/Dropbox/VenomCap_test_data/venom_loci_updated_Mar12_cdhit95_duplicate_exons_renamed_Feb2023_FINAL.fa"
      # output.directory = "/Volumes/LaCie/VenomCap/data-analysis/genetic-distance"
      # mafft.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

      # quiet = TRUE
      # overwrite = FALSE
      # threads = 6
      # memory = 6

      ####################################################################
      ##### Required program path check
      ####################################################################

      require(foreach)

      # Same adds to bbmap path
      if (is.null(mafft.path) == FALSE) {
            b.string <- unlist(strsplit(mafft.path, ""))
            if (b.string[length(b.string)] != "/") {
                  mafft.path <- paste0(append(b.string, "/"), collapse = "")
            } # end if
      } else {
            mafft.path <- ""
      }

      ####################################################################
      ##### Input path check
      ####################################################################

      if (is.null(alignment.directory) == TRUE) {
            stop("Please provide input directory.")
      }

      # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
      if (file.exists(output.directory) == F) {
            dir.create(output.directory)
      } else {
            if (overwrite == TRUE) {
                  system(paste0("rm -r ", output.directory))
                  dir.create(output.directory)
            }
      } # end else

      ####################################################################
      ##### Raw data collection for each alignment and sample
      ####################################################################

      # Gather alignments
      alignment.files <- list.files(alignment.directory, full.names = T)
      alignment.files <- alignment.files[grep(".phy$", alignment.files)]
      target.markers <- Biostrings::readDNAStringSet(target.file)
      target.markers <- target.markers[duplicated(names(target.markers)) != T]

      marker.names <- gsub(".*\\/", "", alignment.files)
      marker.names <- gsub(".phy$", "", marker.names)

      target.markers <- target.markers[names(target.markers) %in% marker.names]
      marker.names <- marker.names[marker.names %in% names(target.markers)]

      # Finds taxa
      taxa.temp <- c()
      for (i in 1:length(alignment.files)) {
            align <- Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = alignment.files[i], format = "phylip"))
            taxa.temp <- unique(append(taxa.temp, names(align)))
      }

      sample.names <- unique(taxa.temp)

      #Sets up multiprocessing
      cl = parallel::makeCluster(threads, outfile = "")
      doParallel::registerDoParallel(cl)
      mem.cl = floor(memory/threads)

      #Loops through each locus and does operations on them
      #all.data = c()
      all.data = foreach::foreach(i=1:length(marker.names),  .combine = rbind, .packages = c("PhyloProcessR", "foreach", "Biostrings","data.table", "ape", "stringr")) %dopar% {
      #Loops through each alignment
      #for (i in 1:length(marker.names)){

      #START HERE
      align.file = alignment.files[grep(paste0(marker.names[i], ".phy$"), alignment.files)]
      align = Biostrings::DNAStringSet(Biostrings::readAAMultipleAlignment(file = align.file, format = "phylip"))

      temp.target = target.markers[names(target.markers) == marker.names[i],]
      temp.target = temp.target[1]
      names(temp.target) = "Reference_Marker"

      #Aligns and then reverses back to correction orientation
      alignment <- PhyloProcessR::runMafft(
            sequence.data = align,
            add.contigs = temp.target,
            algorithm = "add",
            adjust.direction = TRUE,
            threads = 1,
            cleanup.files = T,
            quiet = quiet,
            mafft.path = mafft.path
      )

      #Aligns and then reverses back to correction orientation
      reversed = names(alignment)[grep(pattern = "^_R_", names(alignment))]
      if (length(reversed[grep(pattern = "Reference_Marker", reversed)]) == 1){ alignment = Biostrings::reverseComplement(alignment) }
      names(alignment) = gsub(pattern = "^_R_", replacement = "", x = names(alignment))

      # #Removes the edge gaps
      ref.aligned <- as.character(alignment["Reference_Marker"])
      not.gaps <- stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][, 1]
      ref.start <- min(not.gaps)
      ref.finish <- max(not.gaps)

      # Finds weird gaps to fix
      temp.gaps <- as.numeric(1)
      for (k in 1:length(not.gaps) - 1) {
            temp.gaps <- append(temp.gaps, not.gaps[k + 1] - not.gaps[k])
      }
      temp.gaps <- temp.gaps - 1
      names(temp.gaps) <- not.gaps
      bad.gaps <- which(temp.gaps >= 30)
      front.gaps <- bad.gaps[bad.gaps <= length(not.gaps) * 0.10]
      end.gaps <- bad.gaps[bad.gaps >= length(not.gaps) * 0.90]

      # Fix big gaps if there are any
      if (length(front.gaps) != 0) {
            temp.change <- (max(as.numeric(names(front.gaps)) - ref.start)) - (max(front.gaps) - 1)
            ref.start <- ref.start + temp.change
      } # end gap if

      # Fix big gaps if there are any
      if (length(end.gaps) != 0) {
            add.bp <- length(temp.gaps) - min(end.gaps)
            # add.bp<-(ref.finish-min(as.numeric(names(end.gaps))))
            min.gaps <- temp.gaps[min(end.gaps)]
            temp.change <- as.numeric(names(min.gaps)) - as.numeric(min.gaps)
            ref.finish <- temp.change + add.bp
      } # end gap if

      target.region <- Biostrings::subseq(alignment, ref.start, ref.finish)

      # Removes the edge gaps
      ref.aligned <- as.character(target.region["Reference_Marker"])
      not.gaps <- stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][, 1]
      ref.start <- min(not.gaps)
      ref.finish <- max(not.gaps)
      trim.align <- Biostrings::subseq(target.region, ref.start, ref.finish)

      trim.cols <- PhyloProcessR::trimAlignmentColumns(
            alignment = trim.align,
            min.gap.percent = 50
      )

      #Gets lengths
      write.temp = strsplit(as.character(trim.cols), "")
      aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
      char.align = as.character(as.list(aligned.set))

      #ref.seq = char.align[names(char.align) == "Reference_Marker"]
      #char.align = char.align[names(char.align) != "Reference_Marker"]
      #temp.seq = char.align[1]

      dna.dists = ape::dist.dna(x = aligned.set, model = "raw", as.matrix = TRUE, pairwise.deletion = TRUE)
      dna.frame = data.frame(sample = colnames(dna.dists), dist = dna.dists[rownames(dna.dists) == "Reference_Marker"])
      dna.frame = dna.frame[dna.frame$sample != "Reference_Marker", ]
      
      #Sets up data to collect
      header.all = c("sample" , "target_marker", "target_distance")
      sample.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.all)))
      data.table::setnames(sample.data, header.all)
      sample.data[, sample:=as.character(sample.names)]
      sample.data[, target_marker:=as.character(target_marker)]
      sample.data[, target_marker:=as.character(marker.names[i])]

      data.table::set(sample.data, i = match(dna.frame$sample, sample.data$sample), j = match("target_distance", header.all), value = dna.frame$dist)
      
      sample.data[sample.data$target_distance == 0]$target_distance = NA

      #Writes sample table
      #all.data = rbind(all.data, sample.data)
      print(data.frame(sample.data))

      }#end i loop

      parallel::stopCluster(cl)

      #Writes final table
      write.table(all.data, file = paste0(output.directory, "/target-distance_raw-data.txt"), sep = "\t", row.names = F)

      print(paste0("Saved raw data at ", output.directory, "/target-distance_raw-data.txt"))

      ####################################################################
      ##### Summarizes the raw data previously collected into a summary table
      ####################################################################

      #Get stats
      temp.min = unlist(lapply( split(all.data$target_distance, all.data$sample), min, na.rm = TRUE))
      temp.max = unlist(lapply(split(all.data$target_distance, all.data$sample), max, na.rm = TRUE))
      temp.mean = unlist(lapply(split(all.data$target_distance, all.data$sample), mean, na.rm = TRUE))
      temp.median = unlist(lapply(split(all.data$target_distance, all.data$sample), median, na.rm = TRUE))
      temp.sd = unlist(lapply(split(all.data$target_distance, all.data$sample), sd, na.rm = TRUE))

      #Sets up data to collect
      header.summ = c("sample", "mean_distance", "median_distance",
                  "min_distance", "max_distance", "sd_distance")
      summary.data = data.table::data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.summ)))
      data.table::setnames(summary.data, header.summ)
      summary.data[, sample:=as.character(sample)]
      summary.data[, sample:=as.character(sample.names)]

      data.table::set(summary.data, i = match(names(temp.min), summary.data$sample), j = match("min_distance", header.summ), value = temp.min )
      data.table::set(summary.data, i = match(names(temp.max), summary.data$sample), j = match("max_distance", header.summ), value = temp.max )
      data.table::set(summary.data, i = match(names(temp.mean), summary.data$sample), j = match("mean_distance", header.summ), value = temp.mean )
      data.table::set(summary.data, i = match(names(temp.median), summary.data$sample), j = match("median_distance", header.summ), value = temp.median )
      data.table::set(summary.data, i = match(names(temp.sd), summary.data$sample), j = match("sd_distance", header.summ), value = temp.sd )

      write.table(summary.data, file = paste0(output.directory, "/sample-distance_summary.txt"), sep = "\t", row.names = F)

      print(paste0("Saved raw data at ", output.directory, "/sample-distance_summary.txt"))

}#end function


# ###############################################################################
# ###############################################################################
# ############### Figure 3 Missing Data  w/ markers in common  #################
# ###############################################################################
# ###############################################################################

# #START
# sample.sheet = read_excel("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_April182019.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]

# #Sample directories
# align.dir = "/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Post_Processing_New/Final_Alignments"
# data.scale = list.dirs(align.dir)
# align.folders = data.scale[grep("all-markers_trimmed", data.scale)]
# temp.loci = unique(list.files(align.folders))
# temp.loci = temp.loci[grep("phy", temp.loci)]
# scale.summary<-data.frame()

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")

# ### Hyloidea  
# hyl.loc = fread("Hyloidea_species_sensitivity_raw.txt")
# hyl.dist = fread("Hyloidea_species_distances_raw.txt")

# ### Ranoidea  
# ran.loc = fread("Ranoidea_species_sensitivity_raw.txt")
# ran.dist = fread("Ranoidea_species_distances_raw.txt")

# #Filters
# #filt.ran.loc = ran.loc[ran.loc$Locus %in% hyl.loc$Locus,]
# temp = which(colnames(ran.dist) %in% colnames(hyl.dist))
# filt.ran.dist = ran.dist[,..temp]
# filt.ran.loc = ran.loc[,..temp]
# loci.keep = colnames(filt.ran.loc)

# #Mean of genetic distance from marker per sample
# dist.mean = rowMeans(filt.ran.dist[,4:ncol(filt.ran.dist)], na.rm = TRUE)
# names(dist.mean) = as.character(data.frame(filt.ran.dist[,3])[,1])
# temp.data1 = data.frame(Sample = names(dist.mean), Distance = dist.mean)

# #Mean of genetic distance from marker per sample
# loc.mean = rowMeans(filt.ran.loc[,4:ncol(filt.ran.loc)], na.rm = TRUE)
# names(loc.mean) = as.character(data.frame(filt.ran.loc[,3])[,1])
# temp.data2 = data.frame(Sample = names(loc.mean), MissingBP = 1-loc.mean)

# ran.data = merge(temp.data1, temp.data2, by = "Sample")

# #Filters
# #filt.hyl.loc = hyl.loc[hyl.loc$Locus %in% ran.loc$Locus,]
# temp = which(colnames(hyl.dist) %in% colnames(ran.dist))
# filt.hyl.dist = hyl.dist[,..temp]
# filt.hyl.loc = hyl.loc[,..temp]

# #Mean of genetic distance from marker per sample
# dist.mean = rowMeans(filt.hyl.dist[,4:ncol(filt.hyl.dist)], na.rm = TRUE)
# names(dist.mean) = as.character(data.frame(filt.hyl.dist[,3])[,1])
# temp.data1 = data.frame(Sample = names(dist.mean), Distance = dist.mean)

# #Mean of genetic distance from marker per sample
# loc.mean = rowMeans(filt.hyl.loc[,4:ncol(filt.hyl.loc)], na.rm = TRUE)
# names(loc.mean) = as.character(data.frame(filt.hyl.loc[,3])[,1])
# temp.data2 = data.frame(Sample = names(loc.mean), MissingBP = 1-loc.mean)

# hyl.data = merge(temp.data1, temp.data2, by = "Sample")

# ###################################################
# #For Ranoidea, count samples in each alignment
# ###################################################

# data.combine = c()
# setwd(align.folders[5])
# locus.names = list.files(align.folders[5])
# locus.names = locus.names[gsub(".phy$", "", locus.names) %in% loci.keep]

# sample.names = sample.sheet[sample.sheet$Ranoidea == 1,]$Sample
# comp.matrix = data.frame(Sample = sample.names)

# #Sets up database
# header.data<-c("Sample", locus.names)   
# collect.data<-data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
# setnames(collect.data, header.data)
# collect.data[, Sample:=as.character(Sample)]

# scale.summary = c()
# for (i in 1:length(locus.names)){
  
#   #Read in alignment
#   align = readAAMultipleAlignment(file = locus.names[i], format = "phylip")
  
#   #Use taxa remove
#   tax.names = rownames(align)
#   missing.names = sample.names[!sample.names %in% tax.names]
  
#   set(collect.data, i = match(tax.names, sample.names), j = match(locus.names[i], header.data), 
#       value = as.numeric(1))
  
# } #end i loop

# missing.data = data.frame(Sample = sample.names, Missing = 1-rowSums(collect.data[,2:ncol(collect.data)])/max(rowSums(collect.data[,2:ncol(collect.data)])))
# gen.data = merge(x = ran.data, y = missing.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distance~gen.data$Missing)
# temp.summary<-data.frame(Scale = "Ranoidea Missing Marker", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Missing, y = Distance)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Ranoidea") + xlab("Genetic distance (%)") + ylab("Missing data (% markers)") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p1

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# ggsave("../R_Plots/5_MissingMarker-Dist_Ranoidea.pdf", width = 8, height = 8, useDingbats = F)

# save.data = cbind(Probeset = "Ranoidea", gen.data)
# data.combine = rbind(data.combine, save.data)

# ### DOES relationship for missing bp data

# model<-lm(gen.data$Distance~gen.data$MissingBP)
# temp.summary<-data.frame(Scale = "Ranoidea MissingBP", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)


# ###################################################
# #For Hyloidea, count samples in each alignment
# ###################################################

# setwd(align.folders[3])
# locus.names = list.files(align.folders[3])
# locus.names = locus.names[gsub(".phy$", "", locus.names) %in% loci.keep]

# sample.names = sample.sheet[sample.sheet$Hyloidea == 1,]$Sample
# comp.matrix = data.frame(Sample = sample.names)

# missing.data = data.frame(Sample = sample.names, Missing = 1-rowSums(collect.data[,2:ncol(collect.data)])/max(rowSums(collect.data[,2:ncol(collect.data)])))
# gen.data = merge(x = hyl.data, y = missing.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distance~gen.data$Missing)
# temp.summary<-data.frame(Scale = "Hyloidea Missing Marker", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Missing, y = Distance)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Hyloidea") + xlab("Genetic distance (%)") + ylab("Missing data (% markers)") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p1

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# ggsave("../R_Plots/5_MissingMarker-Dist_Hyloidea.pdf", width = 8, height = 8, useDingbats = F)

# save.data = cbind(Probeset = "Hyloidea", gen.data)
# data.combine = rbind(data.combine, save.data)

# ### DOES relationship for missing bp data

# model<-lm(gen.data$Distance~gen.data$MissingBP)
# temp.summary<-data.frame(Scale = "Hyloidea MissingBP", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)


# ##### GATHER STATS?

# #Quick stats

# ran.data = data.combine[data.combine$ProbeSet %in% "Ranoidea",]
# mean(ran.data$Missing)
# sd(ran.data$Missing)
# min(ran.data$Missing)
# max(ran.data$Missing)

# ran.data = data.combine[data.combine$ProbeSet %in% "Hyloidea",]
# mean(ran.data$Missing)
# sd(ran.data$Missing)
# min(ran.data$Missing)
# max(ran.data$Missing)

# #without crazy samps
# ran.data = ran.data[-21,]
# ran.data = ran.data[-14,]
# ran.data = ran.data[-5,]

# model<-lm(ran.data$Distances~ran.data$Missing)
# temp.summary<-data.frame(Scale = "Hyloidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)


# ran.data = data.combine[data.combine$ProbeSet %in% "Reduced",]
# mean(ran.data$Missing)
# sd(ran.data$Missing)
# min(ran.data$Missing)
# max(ran.data$Missing)


# ###############################################################################
# ###############################################################################
# ########### FIgure: Sample sensivitity with morkers in common #################
# ###############################################################################
# ###############################################################################

# #Sample sensitivity is the amount each target locus is covered by contigs 
# #### *** TEST FOR RELATIONSHIP between sample sensivity and genetic distance

# #START
# sample.sheet<-read_excel("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_April182019.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")

# #############################
# #### A. Sensitivity
# #############################

# scale.summary<-data.frame()

# ### Hyloidea  
# hyl.loc = fread("Hyloidea_species_sensitivity_raw.txt")
# hyl.dist = fread("Hyloidea_species_distances_raw.txt")

# ### Ranoidea  
# ran.loc = fread("Ranoidea_species_sensitivity_raw.txt")
# ran.dist = fread("Ranoidea_species_distances_raw.txt")

# #Filters
# #filt.ran.loc = ran.loc[ran.loc$Locus %in% hyl.loc$Locus,]
# temp = which(colnames(ran.dist) %in% colnames(hyl.dist))
# filt.ran.dist = ran.dist[,..temp]
# filt.ran.loc = ran.loc[,..temp]

# #Mean of genetic distance from marker per sample
# dist.mean = rowMeans(filt.ran.dist[,4:ncol(filt.ran.dist)], na.rm = TRUE)
# names(dist.mean) = as.character(data.frame(filt.ran.dist[,3])[,1])
# temp.data1 = data.frame(Sample = names(dist.mean), Distance = dist.mean)

# #Mean of genetic distance from marker per sample
# loc.mean = rowMeans(filt.ran.loc[,4:ncol(filt.ran.loc)], na.rm = TRUE)
# names(loc.mean) = as.character(data.frame(filt.ran.loc[,3])[,1])
# temp.data2 = data.frame(Sample = names(loc.mean), Locus.Value = loc.mean)

# ran.data = merge(temp.data1, temp.data2, by = "Sample")

# #Filters
# #filt.hyl.loc = hyl.loc[hyl.loc$Locus %in% ran.loc$Locus,]
# temp = which(colnames(hyl.dist) %in% colnames(ran.dist))
# filt.hyl.dist = hyl.dist[,..temp]
# filt.hyl.loc = hyl.loc[,..temp]

# #Mean of genetic distance from marker per sample
# dist.mean = rowMeans(filt.hyl.dist[,4:ncol(filt.hyl.dist)], na.rm = TRUE)
# names(dist.mean) = as.character(data.frame(filt.hyl.dist[,3])[,1])
# temp.data1 = data.frame(Sample = names(dist.mean), Distance = dist.mean)

# #Mean of genetic distance from marker per sample
# loc.mean = rowMeans(filt.hyl.loc[,4:ncol(filt.hyl.loc)], na.rm = TRUE)
# names(loc.mean) = as.character(data.frame(filt.hyl.loc[,3])[,1])
# temp.data2 = data.frame(Sample = names(loc.mean), Locus.Value = loc.mean)

# hyl.data = merge(temp.data1, temp.data2, by = "Sample")
# #hyl.data = hyl.data[hyl.data$Sample != "Bolitaglossa_palmata_LAC_1154",]



# #### DOES THE REGRESSIONS

# model<-lm(ran.data$Locus.Value~ran.data$Distance)
# temp.summary<-data.frame(Scale = "Species", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(ran.data, aes(x = Distance, y = Locus.Value)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Ranoidea Probe Set") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/1_Ranoidea-ProbeSet_sens-dist.pdf", width = 6, height = 6, useDingbats = F)



# model<-lm(hyl.data$Locus.Value~hyl.data$Distance)
# temp.summary<-data.frame(Scale = "Species", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(hyl.data, aes(x = Distance, y = Locus.Value)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Hyloidea Probe Set") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/2_Hyloidea-PRobeSet_sens-dist.pdf", width = 6, height = 6, useDingbats = F)


# ###############################################################################
# ###############################################################################
# ########### FIgure: Sample specificity with morkers in common #################
# ###############################################################################
# ###############################################################################

# #START
# sample.sheet<-read_excel("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_April182019.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")

# #############################
# #### B. Specificity
# #############################

# scale.summary<-data.frame()

# ### Hyloidea  
# all.spec = fread("../All_specificity_summary_raw_probeset_sept9.txt")
# red.spec = fread("../All_specificity_summary_raw_probeset-reduced_sept11.txt")

# ### SUBTRACT FROM REDUCED SPEC
# temp.spec = merge(all.spec, red.spec, by = "Sample")
# new.um = temp.spec$Unmapped.y - (temp.spec$Mapped.x - temp.spec$Mapped.y)

# new.spec = data.frame(ProbeSet = temp.spec$ProbeSet.x, Sample = temp.spec$Sample, 
#                       Mapped = temp.spec$Mapped.y, Unmapped = new.um, Total = temp.spec$Mapped.y + new.um)

# new.spec$PerMap = new.spec$Mapped/new.spec$Total

# #Genetic distances
# hyl.dist = fread("Hyloidea_species_distances_raw.txt")
# ### Ranoidea  
# ran.dist = fread("Ranoidea_species_distances_raw.txt")

# ran.spec = new.spec[new.spec$ProbeSet == "Ranoidea-V1",]
# dist.mean = rowMeans(ran.dist[,4:ncol(ran.dist)], na.rm = TRUE)
# names(dist.mean) = as.character(data.frame(ran.dist[,3])[,1])
# ran.temp = data.frame(Sample = names(dist.mean), Distance = dist.mean)
# ran.data = merge(ran.temp, ran.spec, by = "Sample")

# #Hyloidea
# hyl.spec = new.spec[new.spec$ProbeSet == "Hyloidea-V1",]
# dist.mean = rowMeans(hyl.dist[,4:ncol(hyl.dist)], na.rm = TRUE)
# names(dist.mean) = as.character(data.frame(hyl.dist[,3])[,1])
# hyl.temp = data.frame(Sample = names(dist.mean), Distance = dist.mean)
# hyl.data = merge(hyl.temp, hyl.spec, by = "Sample")
# hyl.data = hyl.data[hyl.data$Sample != "Bolitaglossa_palmata_LAC_1154",]

# # 
# # #Filters
# # #filt.ran.loc = ran.loc[ran.loc$Locus %in% hyl.loc$Locus,]
# # temp = which(colnames(ran.dist) %in% colnames(hyl.dist))
# # filt.ran.dist = ran.dist[,..temp]
# # filt.ran.loc = ran.loc[,..temp]
# # 
# # #Mean of genetic distance from marker per sample
# # dist.mean = rowMeans(filt.ran.dist[,4:ncol(filt.ran.dist)], na.rm = TRUE)
# # names(dist.mean) = as.character(data.frame(filt.ran.dist[,3])[,1])
# # temp.data1 = data.frame(Sample = names(dist.mean), Distance = dist.mean)
# # 
# # #Mean of genetic distance from marker per sample
# # loc.mean = rowMeans(filt.ran.loc[,4:ncol(filt.ran.loc)], na.rm = TRUE)
# # names(loc.mean) = as.character(data.frame(filt.ran.loc[,3])[,1])
# # temp.data2 = data.frame(Sample = names(loc.mean), Locus.Value = loc.mean)
# # 
# # ran.data = merge(temp.data1, temp.data2, by = "Sample")
# # 
# # #Filters
# # #filt.hyl.loc = hyl.loc[hyl.loc$Locus %in% ran.loc$Locus,]
# # temp = which(colnames(hyl.dist) %in% colnames(ran.dist))
# # filt.hyl.dist = hyl.dist[,..temp]
# # filt.hyl.loc = hyl.loc[,..temp]
# # 
# # #Mean of genetic distance from marker per sample
# # dist.mean = rowMeans(filt.hyl.dist[,4:ncol(filt.hyl.dist)], na.rm = TRUE)
# # names(dist.mean) = as.character(data.frame(filt.hyl.dist[,3])[,1])
# # temp.data1 = data.frame(Sample = names(dist.mean), Distance = dist.mean)
# # 
# # #Mean of genetic distance from marker per sample
# # loc.mean = rowMeans(filt.hyl.loc[,4:ncol(filt.hyl.loc)], na.rm = TRUE)
# # names(loc.mean) = as.character(data.frame(filt.hyl.loc[,3])[,1])
# # temp.data2 = data.frame(Sample = names(loc.mean), Locus.Value = loc.mean)
# # 
# # hyl.data = merge(temp.data1, temp.data2, by = "Sample")

# #### DOES THE REGRESSIONS

# model<-lm(ran.data$PerMap~ran.data$Distance)
# temp.summary<-data.frame(Scale = "Species", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(ran.data, aes(x = Distance, y = PerMap)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Ranoidea Probe Set") + ylab("Sample Specificity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)
# p1

# ggsave("../R_Plots/1_Ranoidea-ProbeSet_spec-dist.pdf", width = 6, height = 6, useDingbats = F)



# model<-lm(hyl.data$PerMap~hyl.data$Distance)
# temp.summary<-data.frame(Scale = "Species", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(hyl.data, aes(x = Distance, y = PerMap)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Hyloidea Probe Set") + ylab("Sample Specificity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p1

# ggsave("../R_Plots/2_Hyloidea-ProbeSet_spec-dist.pdf", width = 6, height = 6, useDingbats = F)




# ###############################################################################
# ###############################################################################
# ############### Figure 6: Sample  sensitivity                 #################
# ###############################################################################
# ###############################################################################

# #Panel row 1: Ranoidea, Hyloidea, Reduced Ranoidea 


# #Done Oct 9 2018

# #Get the sample data
# sample.sheet<-read_excel("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_April182019.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]
# scale.names<-colnames(sample.sheet[9:15])

# #Load in data
# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# spp.data<-fread("../All_specificity_summary_raw.txt")
# scale.summary = c()

# ######## Sensitivity
# #######################

# ### Ranoiedea LEVEL 
# sens.data<-fread("Ranoidea_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$Sensitivity)
# temp.summary<-data.frame(Scale = "Ranoidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Sensitivity, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Ranoidea") + xlab("Sample Sensitivity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p1
# ggsave("../R_Plots/6_Sensitivity-Dist_Ranoidea.pdf", width = 8, height = 8)


# ### Hyloidea LEVEL 
# sens.data<-fread("Hyloidea_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$Sensitivity)
# temp.summary<-data.frame(Scale = "Hyloidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Sensitivity, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Hyloidea") + xlab("Sample Sensitivity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p2<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p2
# ggsave("../R_Plots/6_Sensitivity-Dist_Hyloidea.pdf", width = 8, height = 8)


# ######## Specificity
# #######################

# ### Ranoiedea LEVEL 
# sens.data<-fread("Ranoidea_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Ranoidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Ranoidea") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p3<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p3
# ggsave("../R_Plots/6_Spec-Dist_Ranoidea.pdf", width = 8, height = 8)


# ### Hyloidea LEVEL 
# sens.data<-fread("Hyloidea_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Hyloidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Hyloidea") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p4<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p4
# ggsave("../R_Plots/6_Spec-Dist_Hyloidea.pdf", width = 8, height = 8)


# ######## Specificity
# #######################

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# ggsave("R_Plots/6_BoxPlot_Sens_Scale.pdf", width = 9, height = 6)





# #Save combined-all plot
# pdf("../R_Plots/3_ProbeSets_Combined.pdf", width = 10, height = 10)
# multiplot(p1, p2, p3, p4, cols=2)
# dev.off()



# ###############################################################################
# ###############################################################################
# ############### Supplemental: Sample sensivitiy  SCALE        #################
# ###############################################################################
# ###############################################################################

# #Sample sensitivity is the amount each target locus is covered by contigs 
# #### *** TEST FOR RELATIONSHIP between sample sensivity and genetic distance

# #START
# sample.sheet<-read.xls("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_Aug162018.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")

# #############################
# #### A. Phylogenetic Scale
# #############################

# scale.summary<-data.frame()

# ### SPECIES LEVEL 
# gen.data<-fread("Species_species_sens-dist_summary.txt")
# model<-lm(gen.data$Sensitivity~gen.data$Distances)
# temp.summary<-data.frame(Scale = "Species", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Distances, y = Sensitivity)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Species (Cornufer sp.)") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/1_Species_sens-dist.pdf", width = 6, height = 6)


# ### GENUS LEVEL 
# gen.data<-fread("Cornufer_species_sens-dist_summary.txt")
# model<-lm(gen.data$Sensitivity~gen.data$Distances)
# temp.summary<-data.frame(Scale = "Genus", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Distances, y = Sensitivity)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Species (Cornufer sp.)") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p2<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/1_Genus_sens-dist.pdf", width = 6, height = 6)

# ### FAMILY LEVEL 
# gen.data<-fread("Mantellidae_species_sens-dist_summary.txt")
# model<-lm(gen.data$Sensitivity~gen.data$Distances)
# temp.summary<-data.frame(Scale = "Family", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Distances, y = Sensitivity)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Species (Cornufer sp.)") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p3<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/1_Family_sens-dist.pdf", width = 6, height = 6)


# ### RANOIDEA LEVEL 
# gen.data<-fread("Ranoidea_species_sens-dist_summary.txt")
# model<-lm(gen.data$Sensitivity~gen.data$Distances)
# temp.summary<-data.frame(Scale = "Ranoidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Distances, y = Sensitivity)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Species (Cornufer sp.)") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p4<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/1_Ranoidea_sens-dist.pdf", width = 6, height = 6)


# ### HYLOIDEA LEVEL 
# gen.data<-fread("Hyloidea_species_sens-dist_summary.txt")
# model<-lm(gen.data$Sensitivity~gen.data$Distances)
# temp.summary<-data.frame(Scale = "Hyloidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Distances, y = Sensitivity)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Species (Cornufer sp.)") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p5<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/1_Hyloidea_sens-dist.pdf", width = 6, height = 6)


# ### REDUCED LEVEL 
# gen.data<-fread("Reduced_species_sens-dist_summary.txt")
# model<-lm(gen.data$Sensitivity~gen.data$Distances)
# temp.summary<-data.frame(Scale = "Reduced", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Distances, y = Sensitivity)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Species (Cornufer sp.)") + ylab("Sample Sensitivity") + xlab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p6<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/1_Reduced_sens-dist.pdf", width = 6, height = 6)

# #Save combined all plot
# pdf("../R_Plots/1_ALL_sens-dist.pdf", width = 16, height = 20)
# multiplot(p1, p2, p3, p4, p5, p6, cols=2)
# dev.off()

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# write.summary<-cbind(Scale = scale.summary[,1], round(scale.summary[,2:4], 3))
# write.table(write.summary, file = "Scale-dist_stats.txt", row.names = F)


# ###############################################################################
# ###############################################################################
# ############### Supplemental: Sample specificity SCALE        #################
# ###############################################################################
# ###############################################################################

# #Get the sample data
# sample.sheet<-read.xls("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_Aug162018.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]

# #Load in data
# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# spp.data<-fread("../All_specificity_summary_raw.txt")


# ### SPECIES LEVEL 
# sens.data<-fread("Species_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Species", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Species (Cornufer sp.)") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/2_Species_spec-dist.pdf", width = 6, height = 6)


# ### GENUS LEVEL 
# sens.data<-fread("Cornufer_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Genus", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Genus (Cornufer)") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p2<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/2_Genus_spec-dist.pdf", width = 6, height = 6)


# ### Family LEVEL 
# sens.data<-fread("Mantellidae_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Family", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Family (Mantellidae)") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p3<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6)
# #  geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/2_Family_spec-dist.pdf", width = 6, height = 6)


# ### Ranoiedea LEVEL 
# sens.data<-fread("Ranoidea_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Ranoidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Superfamily (Ranoidea)") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p4<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/2_Ranoidea_spec-dist.pdf", width = 6, height = 6)


# ### Hyloidea LEVEL 
# sens.data<-fread("Hyloidea_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Hyloidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Superfamily (Hyloidea)") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p5<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/2_Hyloidea_spec-dist.pdf", width = 6, height = 6)


# ### reduced LEVEL 
# sens.data<-fread("Reduced_species_sens-dist_summary.txt")
# gen.data<-merge(x = sens.data, y = spp.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$PerMap)
# temp.summary<-data.frame(Scale = "Reduced", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = PerMap, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Reduced Ranoidea") + xlab("Sample Specificity") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p6<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# ggsave("../R_Plots/2_Reduced_spec-dist.pdf", width = 6, height = 6)

# #Save combined-all plot
# pdf("../R_Plots/2_ALL_spec-dist.pdf", width = 16, height = 20)
# multiplot(p1, p2, p3, p4, p5, p6, cols=2)
# dev.off()

# write.summary<-cbind(Scale = scale.summary[,1], round(scale.summary[,2:4], 3))
# write.table(write.summary, file = "../Scale-spec_stats.txt", row.names = F)


# #########################################################################################
# ###############################################################################
# ############### Old Versions                                  #################
# ###############################################################################
# #########################################################################################

# ###############################################################################
# ###############################################################################
# ############### Figure 6 -  Specificity summary               #################
# ###############################################################################
# ###############################################################################

# #Get the sample data
# sample.sheet<-read_excel("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_April182019.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]
# scale.names<-colnames(sample.sheet[9:15])

# ###

# #Load in sample file
# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# spp.data<-fread(file = "Species_specificity_summary.txt")

# sample.sheet<-read.xls("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_April182019.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]


# ggplot(data = spp.data, aes(x =ProbeSet, y = PerMap)) + 
#   geom_boxplot(aes(fill = ProbeSet), width = 0.8) + theme_bw()

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# ggsave("R_Plots/5_BoxPlot_Spec_Probeset.pdf", width = 6, height = 4)


# #Loops through each locus and does operations on them
# scale.data<-data.frame()
# for (i in 1:length(scale.names)){
#   #Gets the subdata
#   sub.data<-sample.sheet[sample.sheet[,i+7] == 1,]
#   scale.temp<-spp.data[spp.data$Sample %in% sub.data$Sample,]
#   scale.temp<-cbind(Scale = gsub(".*_", "", scale.names[i]), scale.temp)
#   scale.data<-rbind(scale.data, scale.temp)
# }


# ggplot(data = scale.data, aes(x =Scale, y = PerMap)) + 
#   geom_boxplot(aes(fill = Scale), width = 0.8) + theme_bw()

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# ggsave("R_Plots/5_BoxPlot_Spec_Scale.pdf", width = 9, height = 6)


# ##### Sensivitiy summary
# ###########################

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# a.data<-fread(file = "Species_species_sens-dist_summary.txt")
# b.data<-fread(file = "Reduced_species_sens-dist_summary.txt")
# c.data<-fread(file = "Ranoidea_species_sens-dist_summary.txt")
# d.data<-fread(file = "Hyloidea_species_sens-dist_summary.txt")
# e.data<-fread(file = "Mantellidae_species_sens-dist_summary.txt")
# f.data<-fread(file = "Cornufer_species_sens-dist_summary.txt")
# spp.data<-rbind(a.data, b.data, c.data, d.data, e.data, f.data)

# sample.sheet<-read.xls("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_Aug162018.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]
# scale.names<-unique(sample.sheet$Probe_Set)

# #Loops through each locus and does operations on them
# scale.data<-data.frame()
# for (i in 1:length(scale.names)){
#   #Gets the subdata
#   sub.data<-sample.sheet[sample.sheet$Probe_Set == scale.names[i],]
#   scale.temp<-spp.data[spp.data$Sample %in% sub.data$Sample,]
#   scale.temp<-cbind(ProbeSet = gsub(".*_", "", scale.names[i]), scale.temp)
#   scale.data<-rbind(scale.data, scale.temp)
# }

# ggplot(data = scale.data, aes(x =ProbeSet, y = Sensitivity)) + 
#   geom_boxplot(aes(fill = ProbeSet), width = 0.8) + theme_bw()

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# ggsave("R_Plots/6_BoxPlot_Sens_Probeset.pdf", width = 6, height = 4)


# #Loops through each locus and does operations on them
# scale.names<-colnames(sample.sheet[8:14])
# scale.data<-data.frame()
# for (i in 1:length(scale.names)){
#   #Gets the subdata
#   sub.data<-sample.sheet[sample.sheet[,i+7] == 1,]
#   scale.temp<-spp.data[spp.data$Sample %in% sub.data$Sample,]
#   scale.temp<-cbind(Scale = gsub(".*_", "", scale.names[i]), scale.temp)
#   scale.data<-rbind(scale.data, scale.temp)
# }

# ggplot(data = scale.data, aes(x =Scale, y = Sensitivity)) + 
#   geom_boxplot(aes(fill = Scale), width = 0.8) + theme_bw()

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# ggsave("R_Plots/6_BoxPlot_Sens_Scale.pdf", width = 9, height = 6)




# ###############################################################################
# ###############################################################################
# ############### Figure X Missing Data                         #################
# ###############################################################################
# ###############################################################################

# #1 panel on 1 fig; 4 on another

# # like pis fig w/ prop bp and prop gene overlapping bars

# #1 panel: Alignment missing data across scales AND type of alignment. 

# #4 panels:
# # Genetic distance from before; use %sensivitity subtracted from 1 for proportion. 
# # A. Box plot of gene missing data for the 3 probe sets
# # B. Box plot of base pair missing data for the 3 probe sets
# # C. Hyloidea: missing data vs. genetic distance [done]
# # D. Ranoidea: same  (ONLY gene-wise, bp is reflected in sensitivity) [done]

# ##### Step 1 get sample missing data 

# #START
# sample.sheet = read_excel("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Sample_Rename_Master_Mar132019.xlsx")
# sample.sheet<-sample.sheet[sample.sheet$FrogCap_Paper == 1,]

# #Sample directories
# align.dir = "/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Post_Processing_New/Final_Alignments"
# data.scale = list.dirs(align.dir)
# align.folders = data.scale[grep("all-markers_trimmed", data.scale)]
# temp.loci = unique(list.files(align.folders))
# temp.loci = temp.loci[grep("phy", temp.loci)]

# ###################################################
# #For Ranoidea, count samples in each alignment
# ###################################################

# data.combine = c()
# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# sens.data<-fread("Ranoidea_species_sens-dist_summary.txt")

# setwd(align.folders[5])
# locus.names = list.files(align.folders[5])
# sample.names = sample.sheet[sample.sheet$Ranoidea == 1,]$Sample
# comp.matrix = data.frame(Sample = sample.names)

# #Sets up database
# header.data<-c("Sample", locus.names)   
# collect.data<-data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
# setnames(collect.data, header.data)
# collect.data[, Sample:=as.character(Sample)]

# scale.summary = c()
# for (i in 1:length(locus.names)){
  
#   #Read in alignment
#   align = readAAMultipleAlignment(file = locus.names[i], format = "phylip")
  
#   #Use taxa remove
#   tax.names = rownames(align)
#   missing.names = sample.names[!sample.names %in% tax.names]
  
#   set(collect.data, i = match(tax.names, sample.names), j = match(locus.names[i], header.data), 
#       value = as.numeric(1))
  
# } #end i loop

# missing.data = data.frame(Sample = sample.names, Missing = 1-rowSums(collect.data[,2:ncol(collect.data)])/max(rowSums(collect.data[,2:ncol(collect.data)])))
# gen.data = merge(x = sens.data, y = missing.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$Missing)
# temp.summary<-data.frame(Scale = "Ranoidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Missing, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Ranoidea") + xlab("Missing data") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p1
# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# ggsave("../R_Plots/5_Missing-Dist_Ranoidea.pdf", width = 8, height = 8)

# save.data = cbind("Ranoidea", gen.data)
# data.combine = rbind(data.combine, save.data)

# ###################################################
# #For Hyloidea, count samples in each alignment
# ###################################################

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# sens.data<-fread("Hyloidea_species_sens-dist_summary.txt")

# setwd(align.folders[3])
# locus.names = list.files(align.folders[3])
# sample.names = sample.sheet[sample.sheet$Hyloidea == 1,]$Sample
# comp.matrix = data.frame(Sample = sample.names)

# #Sets up database
# header.data<-c("Sample", locus.names)   
# collect.data<-data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
# setnames(collect.data, header.data)
# collect.data[, Sample:=as.character(Sample)]

# for (i in 1:length(locus.names)){
  
#   #Read in alignment
#   align = readAAMultipleAlignment(file = locus.names[i], format = "phylip")
  
#   #Use taxa remove
#   tax.names = rownames(align)
#   missing.names = sample.names[!sample.names %in% tax.names]
  
#   set(collect.data, i = match(tax.names, sample.names), j = match(locus.names[i], header.data), 
#       value = as.numeric(1))
  
# } #end i loop

# missing.data = data.frame(Sample = sample.names, Missing = 1-rowSums(collect.data[,2:ncol(collect.data)])/max(rowSums(collect.data[,2:ncol(collect.data)])))
# gen.data = merge(x = sens.data, y = missing.data, by.x = "Sample", by.y = "Sample")

# model<-lm(gen.data$Distances~gen.data$Missing)
# temp.summary<-data.frame(Scale = "Hyloidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)

# myplot<-ggplot(gen.data, aes(x = Missing, y = Distances)) + geom_point() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
#         plot.margin = unit(c(1.5, 1, 1, 1), "cm"), #top, right, bottom, left
#         panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
#   ggtitle("Hyloidea") + xlab("Missing data") + ylab("Genetic Distance") +
#   theme(plot.title = element_text(vjust = 10, hjust = 0.5, size = 24)) +
#   theme(axis.title.x = element_text(size = 18, vjust=-0.5),
#         axis.title.y = element_text(size = 18, vjust=1.2))

# p1<-myplot + annotate(geom="text", label=paste0("R2 = ", format(temp.summary$R2, digits = 3), 
#                                                 "; P = ", format(temp.summary$P, digits = 3)),
#                       vjust=1, hjust=0, x = -Inf, y = Inf, size = 6) + geom_smooth(method = "lm", se = T)

# p1
# ggsave("../R_Plots/5_Missing-Dist_Hyloidea.pdf", width = 8, height = 8)

# save.data = cbind("Hyloidea", gen.data)
# data.combine = rbind(data.combine, save.data)

# ##################################################
# #For Reduced-Ranoidea, count samples in each alignment
# ###################################################

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics/Raw_Files")
# sens.data<-fread("Reduced_species_sens-dist_summary.txt")

# setwd(align.folders[6])
# locus.names = list.files(align.folders[6])
# sample.names = sample.sheet[sample.sheet$Reduced == 1,]$Sample
# comp.matrix = data.frame(Sample = sample.names)

# #Sets up database
# header.data<-c("Sample", locus.names)   
# collect.data<-data.table(matrix(as.numeric(0), nrow = length(sample.names), ncol = length(header.data)))
# setnames(collect.data, header.data)
# collect.data[, Sample:=as.character(Sample)]

# for (i in 1:length(locus.names)){
  
#   #Read in alignment
#   align = readAAMultipleAlignment(file = locus.names[i], format = "phylip")
  
#   #Use taxa remove
#   tax.names = rownames(align)
#   missing.names = sample.names[!sample.names %in% tax.names]
  
#   set(collect.data, i = match(tax.names, sample.names), j = match(locus.names[i], header.data), 
#       value = as.numeric(1))
  
# } #end i loop

# missing.data = data.frame(Sample = sample.names, Missing = 1-rowSums(collect.data[,2:ncol(collect.data)])/max(rowSums(collect.data[,2:ncol(collect.data)])))
# gen.data = merge(x = sens.data, y = missing.data, by.x = "Sample", by.y = "Sample")

# save.data = cbind("Reduced", gen.data)
# data.combine = rbind(data.combine, save.data)
# colnames(data.combine) = c("ProbeSet", "Sample", "Sensitivity", "Distances", "Missing")

# #Box plot missing data for each scale
# ggplot(data = data.combine, aes(x =ProbeSet, y = Missing)) + 
#   geom_boxplot(aes(fill = ProbeSet), width = 0.8) + theme_bw()  + ylim(0,1)

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# ggsave("R_Plots/6_BoxPlot_Missing-gene_Probeset.pdf", width = 6, height = 4)

# data.combine$Sensitivity = 1 - data.combine$Sensitivity

# ggplot(data = data.combine, aes(x =ProbeSet, y = Sensitivity)) + 
#   geom_boxplot(aes(fill = ProbeSet), width = 0.8) + theme_bw() + ylim(0,1)

# setwd("/Users/chutter/Dropbox/Research/2_WIP/Anura_Seqcap/Paper_statistics")
# ggsave("R_Plots/6_BoxPlot_Missing-bp_Probeset.pdf", width = 6, height = 4)

# #Quick stats

# ran.data = data.combine[data.combine$ProbeSet %in% "Ranoidea",]
# mean(ran.data$Missing)
# sd(ran.data$Missing)
# min(ran.data$Missing)
# max(ran.data$Missing)

# ran.data = data.combine[data.combine$ProbeSet %in% "Hyloidea",]
# mean(ran.data$Missing)
# sd(ran.data$Missing)
# min(ran.data$Missing)
# max(ran.data$Missing)

# #without crazy samps
# ran.data = ran.data[-21,]
# ran.data = ran.data[-14,]
# ran.data = ran.data[-5,]

# model<-lm(ran.data$Distances~ran.data$Missing)
# temp.summary<-data.frame(Scale = "Hyloidea", R2 = summary(model)$r.squared, Int = summary(model)$coefficients[1,4],
#                          P = summary(model)$coefficients[2,4])
# if (temp.summary$P <= 0.001){ temp.summary$P<-0.001 }
# scale.summary<-rbind(scale.summary, temp.summary)


# ran.data = data.combine[data.combine$ProbeSet %in% "Reduced",]
# mean(ran.data$Missing)
# sd(ran.data$Missing)
# min(ran.data$Missing)
# max(ran.data$Missing)
