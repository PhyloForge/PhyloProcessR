# #REQUIRED PACKAGES
# library(ape)
# library(stringr)
# library(data.table)
# library(GenomicRanges)
# library(Biostrings)
# library(Rsamtools)
# #install.packages("seqinr", dependencies = T, repos = "https://cloud.r-project.org")
# library(seqinr)
#
# options(stringsAsFactors = FALSE)
#
# ##########################################################################################################
# #Parameter setups. Only edit values here.
# ##########################################################################################################
#
# #Set up directories
# threads = "8" #threads, keep quotes
# contig.save = "Rodent"  #This is your save name for the big contig match file
# out.dir = "assembled-contigs" #Output directory name
#
# #CLUSTER directories
# work.dir = "/Volumes/Rodents/Rodents" #Your main project directory
# dir.structure = "analysis" #or "sample" or analysis"
# #proc.dir = "/Volumes/Rodents/Rodents/Processed_Samples" #sample directory
# proc.dir = NA #No sample directory
# contig.dir = "/Volumes/Rodents/Rodents/Test_Contigs"
# prot.file = "/Volumes/Rodents/Selected_Transcripts/Mus_best_prot.fa"
# #prot.file = "/Volumes/Rodents/Selected_Transcripts/protein-mouse_exons-95.fa"
# noncode.file = "/Volumes/Rodents/Selected_Transcripts/Mus_best_noncode.fa"
#
# #The proc dir here is how i have my files organized, with a folder for each sample
# #And all the various steps included in the sample folder. It looks like you have it the other way,
# #Where the steps have their own folders and probably samples inside.
# #I modified the script to take the assembly reads from a folder and set it up like yours w/ NA
# # Processed_Samples
# #   -> Sample_Name
# #     -> raw-reads
# #     -> assembly-reads
# #     -> assembled-contigs
# #     -> etc
#
# ##############################################################################################
# #####################  1.Match loci to contigs                   #############################
# #####################                                            #############################
# ##############################################################################################
#
# #Creates the directory if it doesn't exist
# setwd(work.dir)
# if (dir.structure == "analysis"){
#   if (file.exists(out.dir) == F){ dir.create(out.dir) }
# }#end if
#
# #headers for the blast db
# headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
#            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")
#
# #gets lists of directories and files with sample names
# setwd(contig.dir)
# file.names = list.files(pattern = "", full.names = F, recursive = T)
#
# #Matching and processing for each sample
# for (i in 6:length(file.names)){
#
#   #Sets up working directories for each species
#   sample = gsub(pattern = "_contigs.fa$", replacement = "", x = file.names[i])
#   if (dir.structure == "analysis"){ species.dir = paste0(work.dir, "/", out.dir, "/", sample) }#end if
#   if (dir.structure == "sample"){ species.dir = paste0(proc.dir, "/", sample, "/", out.dir) }#end if
#
#   #Creates species directory if none exists
#   if (file.exists(species.dir) == F){ dir.create(species.dir) }
#   setwd(species.dir)
#   system(paste0("cp ", contig.dir, "/", file.names[i], " ",
#                 species.dir, "/", file.names[i]))
#
#    # DEDUPE almost exact duplicate removal
#   system(paste0("dedupe.sh in=",contig.dir, "/", file.names[i], " ordered=t overwrite=true ",
#                 " out=", sample, "_dedupe.fa", " minidentity=97"), ignore.stderr = T)
#
#   #Make blast database for the probe loci
#    system(paste0("makeblastdb -in ", sample, "_dedupe.fa",
#                  " -parse_seqids -dbtype nucl -out ", sample, "_nucl-blast_db"))
#
#    #Matches samples to proteins
#    system(paste0("tblastn -task tblastn-fast -db ", sample, "_nucl-blast_db -evalue 0.001 -seg no",
#                  " -query ", prot.file, " -out ", sample, "_prot-match.txt",
#                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
#                  " -num_threads ", threads))
#
#   match.data = fread(paste0(sample, "_prot-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
#   setnames(match.data, headers)
#
#   #Matches need to be greater than 12
#   filt.data = match.data[match.data$matches > 12,]
#   #Percent identitiy must match 50% or greater
#   filt.data = filt.data[filt.data$pident >= 50,]
#
#   #Sorting: exon name, contig name, bitscore higher first, evalue
#   setorder(filt.data, qName, tName, -pident, -bitscore, evalue)
#
#   #Make sure the hit is greater than 50% of the reference length
#   filt.data = filt.data[filt.data$matches >= (0.5 * filt.data$qLen),]
#   #Mus: 100,741
#   #Rat: 105,751
#
#   #Mosaic want to find split up targets across a contig or different contigs
#   #find duplicate target matches
#
#   #########################################################################
#   #Part A: Fixes duplicate same contigs matching to the same locus (2 entries)
#   #########################################################################
#
#   contig.names = unique(filt.data[duplicated(filt.data$qName) == T,]$qName)
#
#   #Saves non duplicated data
#   good.data = filt.data[!filt.data$qName %in% contig.names,]
#
#   # Mosaic begin
#   #########################################################################
#   #########################################################################
#   mosaic.count = 0
#   for (j in 1:length(contig.names)) {
#     sub.match = filt.data[filt.data$qName %in% contig.names[j],]
#     sub.match = sub.match[order(sub.match$qStart)]
#
#     if (nrow(sub.match) == 1){ next }
#
#     #Target matches to same contig but broken up
#     if (length(unique(sub.match$tName)) == 1){
#       for (k in 1:(nrow(sub.match)-1)){
#         over.lap = sub.match$qStart[k+1]-sub.match$qEnd[k]
#         if (over.lap > -10){
#           mosaic.count = mosaic.count + 1
#
#         }
#        }#end k
#     }#end if
#
#     #Target matches to same different contigs but broken up
#     if (length(unique(sub.match$tName)) >= 2){
#       for (k in 1:(nrow(sub.match)-1)){
#         over.lap = sub.match$qStart[k+1]-sub.match$qEnd[k]
#         if (over.lap > -10){
#           mosaic.count = mosaic.count + 1
#         }
#       }#end k loop
#     }#end if
#   }#end j
#
#   # Mosaic end
#   #########################################################################
#
#   print(paste0(file.names[i], " ", mosaic.count, " mosaics found!"))
#
# }#end i loop
#
#
#
#
#     #Saves highest bitscore
#     save.match = sub.match[sub.match$bitscore == max(sub.match$bitscore),]
#
#     #Saves longest if equal bitscores
#     save.match = save.match[abs(save.match$qStart-save.match$qEnd) == max(abs(save.match$qStart-save.match$qEnd)),]
#
#     #saves top match here
#     if (nrow(save.match) >= 2){
#       save.match = save.match[1,]
#     }
#
#     new.data = rbind(new.data, save.match)
#
#   } #end j
#
#
#
#
#
#   #Saves final dataset
#   save.data = rbind(good.data, new.data)
#
#   #Reads in contigs
#   contigs = scanFa(FaFile(paste0(sample, "_dedupe.fa")))
#   red.contigs = contigs[names(contigs) %in% save.data$tName]
#   dup.contigs = save.data$tName[duplicated(save.data$tName)]
#   dup.match = save.data[save.data$tName %in% dup.contigs,]
#   dup.data = dup.match[order(dup.match$tName)]
#
#   #Loops through each potential duplicate
#   dup.loci = unique(dup.data$tName)
#   fix.seq = DNAStringSet()
#   for (j in 1:length(dup.loci)){
#     #pulls out data that matches to multiple contigs
#     sub.data = dup.data[dup.data$tName %in% dup.loci[j],]
#     sub.data = sub.data[order(sub.data$tStart)]
#
#     #Fixes direction and adds into data
#     #Finds out if they are overlapping
#     for (k in 1:nrow(sub.data)){
#       new.start = min(sub.data$tStart[k], sub.data$tEnd[k])
#       new.end = max(sub.data$tStart[k], sub.data$tEnd[k])
#       sub.data$tStart[k] = new.start
#       sub.data$tEnd[k] = new.end
#     }#end k loop
#
#     #Saves them if it is split up across the same locus
#     if (length(unique(sub.data$tName)) == 1 && length(unique(sub.data$qName)) == 1){
#       spp.seq = contigs[names(contigs) %in% sub.data$qName]
#       names(spp.seq) = paste0(sub.data$tName[1], "_|_", sample)
#       fix.seq = append(fix.seq, spp.seq)
#       next
#     }
#
#     #Cuts the node apart and saves separately
#     sub.data$tStart = sub.data$tStart-(sub.data$qStart-1)
#     #If it ends up with a negative start
#     sub.data$tStart[sub.data$tStart <= 0] = 1
#     #Fixes ends
#     sub.data$tEnd = sub.data$tEnd+(sub.data$qLen-sub.data$qEnd)
#
#     #Fixes if the contig is smaller than the full target locus
#     sub.data$tEnd[sub.data$tEnd >= sub.data$tLen] = sub.data$tLen[1]
#
#     starts = c()
#     ends = c()
#     starts[1] = 1
#     for (k in 1:(nrow(sub.data)-1)){
#       ends[k] = sub.data$tEnd[k]+floor((sub.data$tStart[k+1]-sub.data$tEnd[k])/2)
#       starts[k+1] = ends[k]+1
#     } #end k loop
#     ends = append(ends, sub.data$tLen[1])
#
#     #Looks for overlapping contigs
#     tmp = ends-starts
#     if(length(tmp[tmp < 0 ]) != 0){
#       starts = sub.data$tStart
#       ends = sub.data$tEnd
#     }
#
#     #Collects new sequence fragments
#     spp.seq = contigs[names(contigs) %in% sub.data$tName]
#     new.seq = DNAStringSet()
#     for (k in 1:length(starts)){ new.seq = append(new.seq, subseq(x = spp.seq, start = starts[k], end = ends[k]) ) }
#
#     #renames and saves
#     names(new.seq) = paste0(sub.data$qName, "_|_", sample)
#     fix.seq = append(fix.seq, new.seq)
#   } #end j loop
#
#   #Writes the base loci
#   base.data = save.data[!save.data$qName %in% gsub("_\\|_.*", "", names(fix.seq)),]
#   base.loci = contigs[names(contigs) %in% base.data$tName]
#   sort.data = base.data[match(names(base.loci), base.data$tName),]
#   #Name and finalize
#   names(base.loci) = paste0(sort.data$qName, "_|_", sample)
#   fin.loci = append(base.loci, fix.seq)
#   fin.loci = fin.loci[width(fin.loci) >= 60]
#
#   #DUPES and numbers don't match up between contigs and table (dupes or not removed?)
#   temp = fin.loci[duplicated(names(fin.loci)) == T]
#   if(length(temp) != 0){ stop("DUPLICATE FOUND") }
#
#   #Finds probes that match to two or more contigs
#   final.loci = as.list(as.character(fin.loci))
#   write.fasta(sequences = final.loci, names = names(final.loci),
#               paste0(sample, "_coding.fa"), nbchar = 1000000, as.string = T)
#
#   merge.contigs = append(merge.contigs, fin.loci)
#
#   system(paste0("rm *nucl-blast_db*"))
#   print(paste0(sample, " probe matching complete. ", length(final.loci), " found!"))
#
# }# end i loop
#
#
# #Saves the coding subset
# setwd(work.dir)
# final.loci = as.list(as.character(merge.contigs))
# write.fasta(sequences = final.loci, names = names(final.loci),
#             paste0("all-sample_coding_contigs.fa"), nbchar = 1000000, as.string = T)
#
#
# ##############################################################################################
# #####################  2. Match noncoding loci to contigs        #############################
# #####################                                            #############################
# ##############################################################################################
#
# #headers for the blast db
# headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
#             "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")
#
# #gets lists of directories and files with sample names
# setwd(contig.dir)
# file.names = list.files(pattern = "", full.names = F, recursive = T)
#
# #Matching and processing for each sample
# merge.contigs = c()
# for (i in 1:length(file.names)){
#
#   #Sets up working directories for each species
#   sample = gsub(pattern = "_contigs.fa$", replacement = "", x = file.names[i])
#   if (dir.structure == "analysis"){ species.dir = paste0(work.dir, "/", out.dir, "/", sample) }#end if
#   if (dir.structure == "sample"){ species.dir = paste0(proc.dir, "/", sample, "/", out.dir) }#end if
#
#   #Creates species directory if none exists
#   if (file.exists(species.dir) == F){ dir.create(species.dir) }
#   setwd(species.dir)
#
#   #Make blast database for the probe loci
#   system(paste0("makeblastdb -in ", sample, "_dedupe.fa",
#                 " -parse_seqids -dbtype nucl -out ", sample, "_nucl-blast_db"))
#
#   #Matches samples to loci
#   system(paste0("blastn -task dc-megablast -db ", sample, "_nucl-blast_db -evalue 0.001",
#                " -query ", noncode.file, " -out ", sample, "_noncode_match.txt",
#                " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
#                " -num_threads ", threads))
#
#   #Need to load in transcriptome for each species and take the matching transcripts to the database
#   match.data = fread(paste0(sample, "_noncode_match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
#   setnames(match.data, headers)
#
#   #Matches need to be greater than 12
#   filt.data = match.data[match.data$matches > 50,]
#   #Percent identitiy must match 50% or greater
#   filt.data = filt.data[filt.data$pident >= 50,]
#
#   #Sorting: exon name, contig name, bitscore higher first, evalue
#   setorder(filt.data, qName, tName, -pident, -bitscore, evalue)
#
#   #Make sure the hit is greater than 50% of the reference length
#   filt.data = filt.data[filt.data$matches >= (0.5 * filt.data$qLen),]
#   #Mus: 100,741
#   #Rat: 105,751
#
#   #########################################################################
#   #Part A: Fixes duplicate same contigs matching to the same locus (2 entries)
#   #########################################################################
#
#   contig.names = unique(filt.data[duplicated(filt.data$qName) == T,]$qName)
#
#   #Saves non duplicated data
#   good.data = filt.data[!filt.data$qName %in% contig.names,]
#
#   new.data = c()
#   save.paralog = c()
#   for (j in 1:length(contig.names)) {
#     sub.match = filt.data[filt.data$qName %in% contig.names[j],]
#     #Skips if 1 row
#     if (nrow(sub.match) == 1){ next }
#
#     #Saves highest bitscore
#     save.match = sub.match[sub.match$bitscore == max(sub.match$bitscore),]
#
#     #Saves longest if equal bitscores
#     save.match = save.match[abs(save.match$qStart-save.match$qEnd) == max(abs(save.match$qStart-save.match$qEnd)),]
#
#     #saves top match here
#     if (nrow(save.match) >= 2){
#       save.match = save.match[1,]
#     }
#
#     new.data = rbind(new.data, save.match)
#
#   } #end j
#
#   #Saves final dataset
#   save.data = rbind(good.data, new.data)
#
#   #Reads in contigs
#   contigs = scanFa(FaFile(paste0(sample, "_dedupe.fa")))
#   red.contigs = contigs[names(contigs) %in% save.data$tName]
#   dup.contigs = save.data$tName[duplicated(save.data$tName)]
#   dup.match = save.data[save.data$tName %in% dup.contigs,]
#   dup.data = dup.match[order(dup.match$tName)]
#
#   #Loops through each potential duplicate
#   dup.loci = unique(dup.data$tName)
#   fix.seq = DNAStringSet()
#   for (j in 1:length(dup.loci)){
#     #pulls out data that matches to multiple contigs
#     sub.data = dup.data[dup.data$tName %in% dup.loci[j],]
#     sub.data = sub.data[order(sub.data$tStart)]
#
#     #Fixes direction and adds into data
#     #Finds out if they are overlapping
#     for (k in 1:nrow(sub.data)){
#       new.start = min(sub.data$tStart[k], sub.data$tEnd[k])
#       new.end = max(sub.data$tStart[k], sub.data$tEnd[k])
#       sub.data$tStart[k] = new.start
#       sub.data$tEnd[k] = new.end
#     }#end k loop
#
#     #Saves them if it is split up across the same locus
#     if (length(unique(sub.data$tName)) == 1 && length(unique(sub.data$qName)) == 1){
#       spp.seq = contigs[names(contigs) %in% sub.data$qName]
#       names(spp.seq) = paste0(sub.data$tName[1], "_|_", sample)
#       fix.seq = append(fix.seq, spp.seq)
#       next
#     }
#
#     #Cuts the node apart and saves separately
#     sub.data$tStart = sub.data$tStart-(sub.data$qStart-1)
#     #If it ends up with a negative start
#     sub.data$tStart[sub.data$tStart <= 0] = 1
#     #Fixes ends
#     sub.data$tEnd = sub.data$tEnd+(sub.data$qLen-sub.data$qEnd)
#
#     #Fixes if the contig is smaller than the full target locus
#     sub.data$tEnd[sub.data$tEnd >= sub.data$tLen] = sub.data$tLen[1]
#
#     starts = c()
#     ends = c()
#     starts[1] = 1
#     for (k in 1:(nrow(sub.data)-1)){
#       ends[k] = sub.data$tEnd[k]+floor((sub.data$tStart[k+1]-sub.data$tEnd[k])/2)
#       starts[k+1] = ends[k]+1
#     } #end k loop
#     ends = append(ends, sub.data$tLen[1])
#
#     #Looks for overlapping contigs
#     tmp = ends-starts
#     if(length(tmp[tmp < 0 ]) != 0){
#       starts = sub.data$tStart
#       ends = sub.data$tEnd
#     }
#
#     #Collects new sequence fragments
#     spp.seq = contigs[names(contigs) %in% sub.data$tName]
#     new.seq = DNAStringSet()
#     for (k in 1:length(starts)){ new.seq = append(new.seq, subseq(x = spp.seq, start = starts[k], end = ends[k]) ) }
#
#     #renames and saves
#     names(new.seq) = paste0(sub.data$qName, "_|_", sample)
#     fix.seq = append(fix.seq, new.seq)
#   } #end j loop
#
#   #Writes the base loci
#   base.data = save.data[!save.data$qName %in% gsub("_\\|_.*", "", names(fix.seq)),]
#   base.loci = contigs[names(contigs) %in% base.data$tName]
#   sort.data = base.data[match(names(base.loci), base.data$tName),]
#   #Name and finalize
#   names(base.loci) = paste0(sort.data$qName, "_|_", sample)
#   fin.loci = append(base.loci, fix.seq)
#   fin.loci = fin.loci[width(fin.loci) >= 60]
#
#   #DUPES and numbers don't match up between contigs and table (dupes or not removed?)
#   temp = fin.loci[duplicated(names(fin.loci)) == T]
#   if(length(temp) != 0){ stop("DUPLICATE FOUND") }
#
#   #Finds probes that match to two or more contigs
#   final.loci = as.list(as.character(fin.loci))
#   write.fasta(sequences = final.loci, names = names(final.loci),
#               paste0(sample, "_noncoding.fa"), nbchar = 1000000, as.string = T)
#
#   merge.contigs = append(merge.contigs, fin.loci)
#
#   system(paste0("rm *nucl-blast_db*"))
#   print(paste0(sample, " probe matching complete. ", length(final.loci), " found!"))
#
# }# end i loop
#
# #Saves the coding subset
# setwd(work.dir)
# final.loci = as.list(as.character(merge.contigs))
# write.fasta(sequences = final.loci, names = names(final.loci),
#             "all-sample_non-coding_contigs.fa", nbchar = 1000000, as.string = T)
#
# ########################################################################
# # STEP 3
# # Output a file of summary stats
# ########################################################################
#
# #gets lists of directories and files with sample names
# setwd(contig.dir)
# file.names = list.files(pattern = "", full.names = F, recursive = T)
#
# #Sets up data summary
# setwd(work.dir)
#
# header.data = c("Sample", "noContigs", "dedupeContigs", "codingTargets", "noncodingTargets",
#                 "totalTargets", "coding_minLen", "coding_maxLen", "coding_meanLen", "coding_medianLen",
#                 "noncoding_minLen", "noncoding_maxLen", "noncoding_meanLen", "noncoding_medianLen")
# samples = gsub("_contigs.fa$", "", file.names)
# samples = gsub(".fa$", "", samples)
#
# prelim.data = data.table(matrix(as.double(0), nrow = length(samples), ncol = length(header.data)))
# setnames(prelim.data, header.data)
# prelim.data[, Sample:=as.character(samples)]
#
# #Cycles through each assembly run and assesses each
# for (i in 1:length(samples)){
#
#   #Sets up working directories for each species
#   if (dir.structure == "analysis"){ species.dir = paste0(work.dir, "/", out.dir, "/", samples[i]) }#end if
#   if (dir.structure == "sample"){ species.dir = paste0(proc.dir, "/", samples[i], "/", out.dir) }#end if
#   setwd(species.dir)
#
#   #Gets length of contigs
#   og.contigs = scanFa(FaFile(paste0(species.dir, "/", samples[i], "_contigs.fa")))
#   dd.contigs = scanFa(FaFile(paste0(species.dir, "/", samples[i], "_dedupe.fa")))
#
#   #Gets the saved matching targets
#   cd.contigs = scanFa(FaFile(paste0(samples[i], "_coding.fa")))
#   nc.contigs = scanFa(FaFile(paste0(samples[i], "_noncoding.fa")))
#
#   set(prelim.data, i = match(samples[i], samples), j = match("noContigs", header.data), value = length(og.contigs) )
#   set(prelim.data, i = match(samples[i], samples), j = match("dedupeContigs", header.data), value = length(dd.contigs) )
#
#   #Gets contig data
#   set(prelim.data, i = match(samples[i], samples), j = match("codingTargets", header.data), value = length(cd.contigs) )
#   set(prelim.data, i = match(samples[i], samples), j = match("noncodingTargets", header.data), value = length(nc.contigs) )
#   set(prelim.data, i = match(samples[i], samples), j = match("totalTargets", header.data), value = length(nc.contigs)+length(cd.contigs) )
#
#   set(prelim.data, i =  match(samples[i], samples), j = match("coding_minLen", header.data), value = min(width(cd.contigs)) )
#   set(prelim.data, i =  match(samples[i], samples), j = match("coding_maxLen", header.data), value = max(width(cd.contigs)) )
#   set(prelim.data, i =  match(samples[i], samples), j = match("coding_meanLen", header.data), value = mean(width(cd.contigs)) )
#   set(prelim.data, i =  match(samples[i], samples), j = match("coding_medianLen", header.data), value = median(width(cd.contigs)) )
#
#   #Noncoding
#   set(prelim.data, i =  match(samples[i], samples), j = match("noncoding_minLen", header.data), value = min(width(nc.contigs)) )
#   set(prelim.data, i =  match(samples[i], samples), j = match("noncoding_maxLen", header.data), value = max(width(nc.contigs)) )
#   set(prelim.data, i =  match(samples[i], samples), j = match("noncoding_meanLen", header.data), value = mean(width(nc.contigs)) )
#   set(prelim.data, i =  match(samples[i], samples), j = match("noncoding_medianLen", header.data), value = median(width(nc.contigs)) )
#
# }#End loop for things
#
# #Saves combined, final dataset
# setwd(work.dir)
# write.csv(prelim.data, file = paste(contig.save, "-sample-assessment.csv", sep = ""))
#
# #### END SCRIPT
# # New files will be located in each species folder within 'Processed_Samples', in the folder assembled_loci
