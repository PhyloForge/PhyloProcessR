#' @title makeExonAlignments
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param alignment.dir path to a folder of sequence alignments in phylip format.
#'
#' @param alignment.format available input alignment formats: fasta or phylip
#'
#' @param output.dir contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.format available output formats: phylip
#'
#' @param HmmCleaner algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param HmmCleaner.path TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param TrimAl if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param TrimAl.path path to a folder of sequence alignments in phylip format.
#'
#' @param trim.external give a save name if you wnat to save the summary to file.
#'
#' @param min.external.percent TRUE to supress mafft screen output
#'
#' @param trim.coverage path to a folder of sequence alignments in phylip format.
#'
#' @param min.coverage.percent contigs are added into existing alignment if algorithm is "add"
#'
#' @param trim.column algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.column.gap.percent TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param alignment.assess if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param min.sample.bp path to a folder of sequence alignments in phylip format.
#'
#' @param min.alignment.length give a save name if you wnat to save the summary to file.
#'
#' @param min.taxa.alignment TRUE to supress mafft screen output
#'
#' @param min.gap.percent if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param overwrite TRUE to supress mafft screen output
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

# makeExonAlignments = function(alignment.directory = NULL,
#                               alignment.format = "phylip",
#                               output.directory = NULL,
#                               output.format = "phylip",
#                               reference.type = c("target", "alignment"),
#                               reference.path = NULL,
#                               target.direction = TRUE,
#                               concatenate.intron.flanks = TRUE,
#                               threads = 1,
#                               memory = 1,
#                               overwrite = FALSE,
#                               resume = TRUE,
#                               mafft.path = NULL) {

#
#   alignment.directory = "alignments/untrimmed_all-markers"
#   alignment.format = "phylip"
#   output.directory = "alignments/untrimmed_introns"
#   output.format = "phylip"
#   reference.type = "target"
#   reference.path = target.file
#   target.direction = TRUE
#   concatenate.intron.flanks = TRUE
#   threads = threads
#   memory = memory
#   overwrite = overwrite
#   resume = resume
#   mafft.path = mafft.path

#
#   #Same adds to bbmap path
#   if (is.null(mafft.path) == FALSE){
#     b.string = unlist(strsplit(mafft.path, ""))
#     if (b.string[length(b.string)] != "/") {
#       mafft.path = paste0(append(b.string, "/"), collapse = "")
#     }#end if
#   } else { mafft.path = "" }
#
#   if (alignment.directory == output.directory){ stop("You should not overwrite the original alignments.") }
#
#   # if (dir.exists(output.dir) == FALSE) { dir.create(output.dir) }
#
#   #So I don't accidentally delete everything while testing resume
#   if (resume == TRUE & overwrite == TRUE){
#     overwrite = FALSE
#     stop("Error: resume = T and overwrite = T, cannot resume if you are going to delete everything!")
#   }
#
#   if (dir.exists(output.directory) == TRUE) {
#     if (overwrite == TRUE){
#       system(paste0("rm -r ", output.directory))
#       dir.create(output.directory)
#     }
#   } else { dir.create(output.directory) }
#
#   #Gathers alignments
#   align.files = list.files(alignment.directory)
#
#   if (reference.type == "target"){
#     target.loci = Biostrings::readDNAStringSet(file = reference.path, format = "fasta")
#   }#end if
#
#   if (reference.type == "alignment"){
#     ref.align = list.files(reference.path)
#   }#end if
#
#   if (length(align.files) == 0) { stop("alignment files could not be found.") }
#
#   #Skips files done already if resume = TRUE
#   if (resume == TRUE){
#     done.files = list.files(output.directory)
#     align.files = align.files[!gsub("\\..*", "", align.files) %in% gsub("\\..*", "", done.files)]
#   }
#
#   if (length(align.files) == 0) { stop("All alignments have already been completed and overwrite = FALSE.") }
#
#   #Sets up multiprocessing
#   cl = parallel::makeCluster(threads, outfile = "")
#   doParallel::registerDoParallel(cl)
#   mem.cl = floor(memory/threads)
#
#   #Loops through each locus and does operations on them
#   foreach(i=1:length(align.files), .packages = c("PhyloCap", "foreach", "Biostrings", "ape", "stringr")) %dopar% {
#   #Loops through each locus and does operations on them
#   #for (i in 1:length(align.files)) {
#     #Load in alignments
#     if (alignment.format == "phylip"){
#       align = Biostrings::readAAMultipleAlignment(file = paste0(alignment.directory, "/", align.files[i]), format = "phylip")
#
#       #  align = Biostrings::readDNAStringSet(file = paste0(alignment.dir, "/", align.files[i]), format = "phylip")
#       #  align = readLines(paste0(alignment.dir, "/", align.files[i]))[-1]
#       #  align = gsub(".*\\ ", "", align)
#       #  char.count = nchar(align)
#
#       align = Biostrings::DNAStringSet(align)
#       save.name = gsub(".phy$", "", align.files[i])
#       save.name = gsub(".phylip$", "", save.name)
#     }#end phylip
#
#     if (alignment.format == "fasta"){
#       align = Biostrings::readDNAStringSet(paste0(alignment.directory, "/", align.files[i]) )
#       save.name = gsub(".fa$", "", align.files[i])
#       save.name = gsub(".fasta$", "", save.name)
#     }#end phylip
#
#     #First steps here?
#     ##### *** START HERE
#
#     trimmed = trim.ends(align, min.n.seq = ceiling(length(align) * 0.51), codon.trim = T)
#     if (length(trimmed) <= min.taxa){ return(NULL) }
#     save.names = names(trimmed)
#
#     #Gets consensus seq for trimming more
#     con.seq = make.consensus(trimmed, method = "majority")
#
#     #Removes the edge gaps
#     ref.aligned = as.character(con.seq)
#     not.gaps = str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
#
#     #Finds weird gaps to fix
#     temp.gaps = as.numeric(1)
#     for (k in 1:length(not.gaps)-1){ temp.gaps = append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
#     temp.gaps = temp.gaps-1
#     names(temp.gaps) = not.gaps
#     gap.spots = temp.gaps[temp.gaps %% 3 != 0]
#
#     del.col<-c()
#     if (length(gap.spots) != 0){
#       #Loop through each potential bad thing and fix
#       for (k in 1:length(gap.spots)){
#         del.col = append(del.col, (as.numeric(names(gap.spots[k]))-gap.spots[k]):(as.numeric(names(gap.spots[k]))-1))
#       }
#     }#end if
#
#     ##############
#     #STEP 3: Fixes large gaps at ends of alignment
#     ##############
#     #Looks for gaps that clearly wrong and not 3 BP
#     new.align = strsplit(as.character(trimmed), "")
#     x = as.matrix(as.DNAbin(new.align))
#
#     rem.n<-c()
#     for (k in 1:ncol(x)){
#       gaps<-table(as.character(x[,k]))
#       per.gaps<-gaps[names(gaps) == "-"]/nrow(x)
#
#       if (length(per.gaps) != 0){
#         #Records column when the gaps exceed this percentage
#         if (per.gaps >= 0.75){ del.col<-append(del.col, k) }
#         #Removes gap columns only consisting of Ns
#         n.gaps<-gaps[names(gaps) != "-"]
#         if (length(n.gaps) == 1){
#           if (names(n.gaps) == "n"){ rem.n<-append(rem.n, k)}
#         }
#       }#end if
#
#     }#end k loop
#
#     #combines columns to be deleted
#     fin.del<-c(rem.n, del.col)
#     if (length(fin.del) != 0){
#       x<-x[,-fin.del] }
#     #Removes bad columsn and coverts alignment back to DNASTringSet
#     char.align<-as.list(data.frame(t(as.character(x))))
#     temp.align<-lapply(char.align, FUN = function(x) paste(x, collapse = ""))
#     trimmed<-DNAStringSet(unlist(temp.align))
#     names(trimmed)<-save.names
#
#     ##############
#     #STEP 4: Gathers table of best and longest stop codon free frames for each seq
#     ##############
#     #Checks to make sure the codon position is correct
#     save.frame = data.frame()
#     save.all = data.frame()
#     for (j in 1:length(trimmed)){
#       #Finds open reading frames
#       temp.codon = find.orf(trimmed[j], mitochondrial = F, min.size = 50 )
#       if(nrow(temp.codon) == 0){
#         samp.frame<-cbind(Sample = names(trimmed[j]), FrameStart = 0, FrameEnd = 0, Size = 0, sppSize = 0, Frame = "0")
#         save.frame<-rbind(save.frame, samp.frame)
#       } else {
#
#         all.frame<-temp.codon[temp.codon$Size >= max(temp.codon$Size) * .70,]
#         big.frame<-temp.codon[temp.codon$Size == max(temp.codon$Size),]
#
#         if (nrow(big.frame) >= 2){
#           #Picks the best from this order of things
#           temp.stop<-big.frame[big.frame$Frame == "F1",]
#           if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R1",] }
#           if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F2",] }
#           if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R2",] }
#           if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "F3",] }
#           if (nrow(temp.stop) == 0){ temp.stop<-big.frame[big.frame$Frame == "R3",] }
#           big.frame<-temp.stop
#         }
#
#         # #Saves teh data
#         samp.frame<-cbind(Sample = names(trimmed[j]), big.frame)
#         temp.size<-unlist(strsplit(as.character(trimmed[j]), ""), use.names = F)
#
#         #Starts from the beginning and end to fill in end gaps
#         sub.size<-0
#         for (q in 1:length(temp.size)){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
#         for (q in length(temp.size):1){ if (temp.size[q] == "-"){ sub.size<-sub.size+1 } else { break } }
#
#         #Saves final data
#         samp.frame<-cbind(samp.frame, sppSize = length(temp.size)-sub.size)
#         save.frame<-rbind(save.frame, samp.frame)
#         all.frame<-cbind(Sample = names(trimmed[j]), all.frame)
#         all.frame<-cbind(all.frame, sppSize = length(temp.size)-sub.size)
#         save.all<-rbind(save.all, all.frame)
#       }#end else
#     }#end j loop
#
#     #Moves on if there are no frames found. Saves to Anon folder?
#     if (unique(save.frame$Frame)[1] == "0"){
#       print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
#       return(NULL)
#     }
#
#     ##############
#     #STEP 5: Uses previous data to find a consistent frame
#     ##############
#     #Looks at the overall data rather than the best indiv data to find a consistent frame
#     temp.all = save.all
#     frame.names = unique(temp.all$Frame)
#     #Goes through the equally good frames and reduces to frames with the same range
#     very.best = data.frame()
#     for (k in 1:length(frame.names)){
#       temp.best<-temp.all[temp.all$Frame == frame.names[k],]
#       starts<-table(temp.best$FrameStart)[table(temp.best$FrameStart) == max(table(temp.best$FrameStart))]
#       ends<-table(temp.best$FrameEnd)[table(temp.best$FrameEnd) == max(table(temp.best$FrameEnd))]
#
#       #Removes duplicates
#       starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
#       ends<-ends[as.numeric(names(ends)) == min(as.numeric(names(ends)))]
#
#       super.best<-temp.best[temp.best$FrameStart == as.numeric(names(starts)),]
#       super.best<-super.best[super.best$FrameEnd == as.numeric(names(ends)),]
#       very.best<-rbind(very.best, super.best)
#     }#end k loop
#
#     #Moves on if there are no frames found. Saves to Anon folder?
#     if (nrow(very.best) == 0){
#       print(paste(locus.names[i], " sucked. No frames found.", sep = ""))
#       return(NULL)
#     }
#
#     ##############
#     #STEP 6: Selects the best frame
#     ##############
#     #Picks out the best frame
#     best.frame = table(very.best$Frame)[table(very.best$Frame) == max(table(very.best$Frame))]
#
#     #If there are multiple good frames pick the biggest
#     if (length(best.frame) != 1){
#       temp.fix<-very.best[very.best$Frame %in% names(best.frame),]
#       bigger<-temp.fix[temp.fix$Size == max(temp.fix$Size),]
#       best.frame<-table(bigger$Frame)[table(bigger$Frame) == max(table(bigger$Frame))]
#     }#end if
#
#     #If they are same size just pick from this order
#     if (length(best.frame) != 1){
#       #Picks the best from this order of things
#       temp.stop<-best.frame[names(best.frame) == "F1"]
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R1"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F2"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R2"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "F3"] }
#       if (length(temp.stop) == 0){ temp.stop<-best.frame[names(best.frame) == "R3"] }
#       best.frame<-temp.stop
#     }
#
#     #Checks the remaining ORF size
#     temp.size<-save.all[save.all$Frame == names(best.frame),]
#     samp.spp<-temp.size[temp.size$sppSize == max(temp.size$sppSize),]
#     samp.seq<-trimmed[names(trimmed) == samp.spp$Sample[1]]
#     samp.size<-nchar(gsub("-", "", as.character(samp.seq)))
#
#     if (mean(temp.size$Size) <= samp.size*.5){ print(paste(locus.names[i], " was small.", sep = "")) }
#
#     if (mean(temp.size$Size) <= samp.size*.25){
#       print(paste(locus.names[i], " sucked. Too few sequence left.", sep = ""))
#       return(NULL)
#     }
#
#     if (best.frame <= length(trimmed) * .5){
#       print(paste(locus.names[i], " sucked. No cosistent frame.", sep = ""))
#       return(NULL)
#     }
#
#     #Reverses if it needs to
#     if (length(grep("R", names(best.frame))) != 0){
#       new.align<-reverseComplement(trimmed)
#     }else { new.align<-trimmed }
#
#     ##############
#     #STEP 7: Gets start and stop coordinates for each sequence and find best alignment
#     ##############
#     #Gets trimming locations
#     frame.ranges<-save.all[save.all$Frame == names(best.frame),]
#
#     #Gets potential starts and ends
#     starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
#     ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
#     starts<-starts[starts == max(starts)]
#     ends<-ends[ends == max(ends)]
#
#     if (length(starts) != 1 || length(ends) != 1){
#       frame.ranges<-save.all[save.all$Frame == names(best.frame),]
#       starts<-table(frame.ranges$FrameStart)[table(frame.ranges$FrameStart) == max(table(frame.ranges$FrameStart))]
#       ends<-table(frame.ranges$FrameEnd)[table(frame.ranges$FrameEnd) == max(table(frame.ranges$FrameEnd))]
#       starts<-starts[starts == max(starts)]
#       ends<-ends[ends == max(ends)]
#     }
#
#     if (length(starts) != 1 || length(ends) != 1){
#       starts<-starts[as.numeric(names(starts)) == min(as.numeric(names(starts)))]
#       ends<-ends[as.numeric(names(ends)) == max(as.numeric(names(ends)))]
#     }
#
#     ###################
#     #STEP 8: Makes sure entire alignment is a multiple of 3
#     ###################
#     anu.start<-as.numeric(names(starts))
#     new.end<-as.numeric(names(ends))
#     new.len<-new.end-(anu.start-1)
#
#     #Gets a new end to keep in multiples of 3 for proteins
#     if (length(new.len[which(new.len %%3==0)]) == 0) {
#       anu.end<-new.end-1
#     } else { anu.end<-new.end }
#
#     new.len<-anu.end-(anu.start-1)
#     if (length(new.len[which(new.len %%3==0)]) == 0) {
#       anu.end<-new.end-2
#     } else { anu.end<-anu.end }
#
#     #Trims sequence with new coords
#     done.seq<-subseq(start = anu.start, end = anu.end, x = new.align)
#
#     ###################
#     #STEP 9: Trim out odd start/end bases
#     ###################
#     codon.seq<-DNAStringSet()
#     for (k in 1:length(done.seq)){
#       ref.aligned<-as.character(done.seq[k])
#
#       #Chcecks at beginning of sequence
#       not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
#
#       if (length(not.gaps) >= min.len){
#
#         #Checks if its odd, delete 1 base
#         if ( (not.gaps[1]-1) %%3 == 2){ substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-" }
#
#         #Deletes 2 bases its off by
#         if ( (not.gaps[1]-1) %%3 == 1){
#           substr(ref.aligned, not.gaps[1], not.gaps[1])<-"-"
#           substr(ref.aligned, not.gaps[1]+1, not.gaps[1]+1)<-"-"
#         }#end if
#
#         #checks for end of sequence
#         not.gaps<-str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
#         #counts characters
#         char.len<-(not.gaps[length(not.gaps)]-not.gaps[1])+1
#         end.pos<-not.gaps[length(not.gaps)]
#         #removes odd characters at ends
#         if ( char.len %%3 == 1){ substr(ref.aligned, end.pos, end.pos)<-"-" }
#
#         if ( char.len %%3 == 2){
#           substr(ref.aligned, end.pos-1, end.pos-1)<-"-"
#           substr(ref.aligned, end.pos, end.pos)<-"-"
#         } #end if
#
#         #Saves final seqs
#         save.seq<-DNAStringSet(ref.aligned)
#         names(save.seq)<-names(done.seq[k])
#         codon.seq<-append(codon.seq, save.seq)
#       }#end if
#
#     } #END K
#
#     ###################
#     #STEP 10: Change stop codons to N
#     ###################
#
#     if (length(codon.seq) != 0){
#       #Finds stop codons to replace
#       n.seq = DNAStringSet(gsub("-", "N", as.character(codon.seq)))
#       stop.seq = DNAStringSet()
#       del.taxa = c()
#       for (k in 1:length(n.seq)){
#         stop.data = find.codon(n.seq[k], genetic.code = "nuclear", reverse = T)
#         stop.data = stop.data[stop.data$frame == "F1",]
#         #If there is a stop codon take bold action
#         if (stop.data$start[1] == 0){
#           stop.seq = append(stop.seq, n.seq[k])
#         } else {
#           #Goes through each codon
#           #  for (y in 1:nrow(stop.data)){
#           #Saves final seqs
#           #    substr(ref.aligned, stop.data$start[y], stop.data$start[y])<-"N"
#           #    substr(ref.aligned, stop.data$start[y]+1, stop.data$start[y]+1)<-"N"
#           #    substr(ref.aligned, stop.data$start[y]+2, stop.data$start[y]+2)<-"N"
#           #  }#end Y LOOP
#           del.taxa = append(del.taxa, names(ref.aligned) )
#           #Saves final data
#           #save.seq = DNAStringSet(ref.aligned)
#           #names(save.seq) = names(n.seq[k])
#           #stop.seq = append(stop.seq, save.seq)
#         }#end if state
#       }# END K loop
#     } else {
#       stop.seq = done.seq
#     }#end else
#
#     if (length(stop.seq) <= min.taxa){ return(NULL) }
#
#     ###################
#     #FINAL STEP: Save everything after some final spp and length filtering
#     ###################
#     #Removes sequences that are less than a certain coverage\
#     t.align<-strsplit(as.character(stop.seq), "")
#     len.loci<-lapply(t.align, function (x) x[x != "N"])
#     spp.len<-unlist(lapply(len.loci, function (x) length(x)))
#     spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]
#     spp.rem<-append(spp.rem, spp.len[spp.len <= as.numeric(min.len)]  )
#
#     if (length(spp.rem) > 0){
#       red.align<-t.align[!names(t.align) %in% unique(names(spp.rem))]
#     } else { red.align<-t.align }
#
#     #Removes if all removed
#     if (length(red.align) == 0){
#       print(paste0("not enough alignment remains for ", locus.names[i]))
#       return(NULL)
#     }
#
#     #removes loci with too few taxa
#     if (length(names(red.align)) <= as.numeric(min.taxa)){
#       print(paste("Too few taxa after trimming."))
#       return(NULL)
#     }
#
#     #writes alignment
#     mat.align<-lapply(red.align, tolower)
#     write.align<-as.matrix(as.DNAbin(mat.align))
#
#     #readies for saving
#     write.phy(write.align, file=paste0("exon-only_trimmed/", locus.names[i]), interleave = F)
#
# }#end function

    ##### ==== OLD
#
#
#     if (reference.type == "target"){
#       #Loads in a pulls out relevant target file
#       target.seq = target.loci[names(target.loci) %in% save.name]
#       names(target.seq) = "Reference_Locus"
#
#     }#end target if
#
#     #If using the alignments
#     if (reference.type == "alignment"){
#
#       if (file.exists(paste0(reference.path, "/", align.files[i])) == FALSE){ return(NULL) }
#
#       ref.align = Biostrings::readAAMultipleAlignment(file = paste0(reference.path, "/", align.files[i]), format = "phylip")
#       ref.align = Biostrings::DNAStringSet(ref.align)
#       #Gets consensus seq for trimming more
#       con.seq = makeConsensus(ref.align, method = "majority")
#       #Removes the edge gaps
#       target.seq = Biostrings::DNAStringSet(gsub("\\+|-", "", as.character(con.seq)))
#       names(target.seq) = "Reference_Locus"
#     }#end alignment if
#
#
#     #Checks for correct target sequence amount
#     if (length(target.seq) == 0){ return(NULL) }
#     if (length(target.seq) >= 2){ return(NULL) }
#     if (Biostrings::width(target.seq) <= 10) { return(NULL) }
#
#     ##############
#     #STEP 2: Runs MAFFT to add
#     ##############
#     #Aligns and then reverses back to correction orientation
#     alignment = runMafft(sequence.data = align,
#                          add.contigs = target.seq,
#                          save.name = paste0(output.directory, "/", align.files[i]),
#                          algorithm = "add",
#                          adjust.direction = TRUE,
#                          threads = 1,
#                          cleanup.files = T,
#                          quiet = TRUE,
#                          mafft.path = mafft.path)
#
#     #Checks for failed mafft run
#     if (length(alignment) == 0){ return(NULL) }
#
#     #Checks if you want to keep to target direction or not
#     if (target.direction == TRUE){
#       #Aligns and then reverses back to correction orientation
#       reversed = names(alignment)[grep(pattern = "_R_", names(alignment))]
#       if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){ alignment = Biostrings::reverseComplement(alignment) }
#       names(alignment) = gsub(pattern = "_R_", replacement = "", x = names(alignment))
#     } else {
#       #Regular fixes
#       names(alignment) = gsub(pattern = "_R_", replacement = "", x = names(alignment))
#     }#end if
#
#     # ##############
#     # #STEP 3: Removes exon from the intron part
#     # ##############
#     # #Removes the edge gaps
#     ref.aligned = as.character(alignment['Reference_Locus'])
#     not.gaps = stringr::str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
#     ref.start = min(not.gaps)
#     ref.finish = max(not.gaps)
#
#     #Finds weird gaps to fix
#     temp.gaps = as.numeric(1)
#     for (k in 1:length(not.gaps)-1){ temp.gaps = append(temp.gaps, not.gaps[k+1]-not.gaps[k]) }
#     temp.gaps = temp.gaps-1
#     names(temp.gaps) = not.gaps
#     bad.gaps = which(temp.gaps >= 30)
#     front.gaps = bad.gaps[bad.gaps <= length(not.gaps) *0.10]
#     end.gaps = bad.gaps[bad.gaps >= length(not.gaps) *0.90]
#
#     #Fix big gaps if there are any
#     if (length(front.gaps) != 0){
#       temp.change = (max(as.numeric(names(front.gaps))-ref.start))-(max(front.gaps)-1)
#       ref.start = ref.start+temp.change
#     }#end gap if
#
#     #Fix big gaps if there are any
#     if (length(end.gaps) != 0){
#       add.bp = length(temp.gaps)-min(end.gaps)
#       #add.bp<-(ref.finish-min(as.numeric(names(end.gaps))))
#       min.gaps = temp.gaps[min(end.gaps)]
#       temp.change = as.numeric(names(min.gaps))-as.numeric(min.gaps)
#       ref.finish = temp.change+add.bp
#     }#end gap if
#
#     #Cuts out the intron pieces
#     intron.left = Biostrings::subseq(alignment, 1, ref.start-1)
#     intron.right = Biostrings::subseq(alignment, ref.finish+1, Biostrings::width(alignment))
#     save.names  = names(alignment)
#
#     #saves intron flanks separately
#     if (concatenate.intron.flanks == FALSE) {
#       #Remove gap only alignments
#       gap.align = strsplit(as.character(intron.left), "")
#       gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
#       gap.rem = gap.count[gap.count <= 10]
#       left.intron = intron.left[!names(intron.left) %in% names(gap.rem)]
#
#       #Remove gap only alignments
#       gap.align = strsplit(as.character(intron.right), "")
#       gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
#       gap.rem = gap.count[gap.count <= 10]
#       right.intron = intron.right[!names(intron.right) %in% names(gap.rem)]
#
#       #Saves them
#       write.temp = strsplit(as.character(left.intron), "")
#       aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
#
#       save.name = gsub(".phy$", "", align.files[i])
#
#       #readies for saving
#       writePhylip(alignment = aligned.set,
#                   file=paste0(output.directory, "/", save.name, "_1.phy"),
#                   interleave = F,
#                   strict = F)
#
#       #Saves them
#       write.temp = strsplit(as.character(right.intron), "")
#       aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
#
#       #readies for saving
#       writePhylip(alignment = aligned.set,
#                   file=paste0(output.directory, "/", save.name, "_2.phy"),
#                   interleave = F,
#                   strict = F)
#
#       #Deletes old files
#       print(paste0(align.files[i], " alignment saved."))
#
#     } else {
#
#       #Merges the alignments
#       intron.align = Biostrings::DNAStringSet(paste0(as.character(intron.left), as.character(intron.right)))
#       names(intron.align) = save.names
#       intron.align = intron.align[names(intron.align) != "Reference_Locus"]
#
#       #Remove gap only alignments
#       gap.align = strsplit(as.character(intron.align), "")
#       gap.count = unlist(lapply(gap.align, function(x) length(x[x != "-"]) ) )
#       gap.rem = gap.count[gap.count <= 10]
#       rem.align = intron.align[!names(intron.align) %in% names(gap.rem)]
#
#       ##############
#       #STEP 5: Cleanup and save
#       ##############
#
#       if (length(rem.align) != 0 && length(Biostrings::width(rem.align)) != 0){
#         #string splitting
#         write.temp = strsplit(as.character(rem.align), "")
#         aligned.set = as.matrix(ape::as.DNAbin(write.temp) )
#
#         #readies for saving
#         writePhylip(alignment = aligned.set,
#                     file=paste0(output.directory, "/", align.files[i]),
#                     interleave = F,
#                     strict = F)
#
#         #Deletes old files
#         print(paste0(align.files[i], " alignment saved."))
#       } else { print(paste0(align.files[i], " alignment not saved. Not enough data."))  }
#     }#end
#
#     rm()
#     gc()
#
#   }#end i loop
#
#   #close multithread
#   parallel::stopCluster(cl)
#
# } #end function
