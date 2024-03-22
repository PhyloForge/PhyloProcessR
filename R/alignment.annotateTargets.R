#' @title annotateTargets
#'
#' @description Function for batch trimming a folder of alignments, with the various trimming functions available to select from
#'
#' @param assembly.directory path to a folder of sequence alignments in phylip format.
#'
#' @param target.file available input alignment formats: fasta or phylip
#'
#' @param alignment.contig.name contigs are added into existing alignment if algorithm is "add"
#'
#' @param output.directory available output formats: phylip
#'
#' @param min.match.percent algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param min.match.length TRUE applies the adjust sequence direction function of MAFFT
#'
#' @param min.match.coverage if a file name is provided, save.name will be used to save aligment to file as a fasta
#'
#' @param threads path to a folder of sequence alignments in phylip format.
#'
#' @param memory give a save name if you wnat to save the summary to file.
#'
#' @param trim.target TRUE to supress mafft screen output
#'
#' @param overwrite path to a folder of sequence alignments in phylip format.
#'
#' @param resume contigs are added into existing alignment if algorithm is "add"
#'
#' @param quiet algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param blast.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param bbmap.path algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
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

annotateTargets = function(assembly.directory = NULL,
                            target.file = NULL,
                            alignment.contig.name = "annotated-contigs-all",
                            output.directory = "annotated-contigs",
                            min.match.percent = 60,
                            min.match.length = 60,
                            min.match.coverage = 50,
                            trim.target = FALSE,
                            threads = 1,
                            memory = 1,
                            blast.path = NULL,
                            cdhit.path = NULL,
                            overwrite = FALSE,
                            quiet = TRUE
                            ) {
#
  #Debug setup
  #Debugging
  # setwd("/Volumes/LaCie/Microhylidae_test/")
  # assembly.directory <- "/Volumes/LaCie/Microhylidae_test/data-analysis/contigs/7_filtered-contigs"
  # target.file = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
  # output.directory = "data-analysis/contigs/8_annotated-contigs"
  # blast.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # cdhit.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  #
  # #
  # # #Main settings
  # threads = 8
  # memory = 40
  # trim.target = FALSE
  # overwrite = FALSE
  # quiet = TRUE
  # #
  # # #tweak settings (make some statements to check these)
  # min.match.percent = 60
  # min.match.length = 70
  # min.match.coverage = 35
  # #

  require(foreach)

  #Add the slash character to path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  # Same adds to bbmap path
  if (is.null(cdhit.path) == FALSE) {
    b.string <- unlist(strsplit(cdhit.path, ""))
    if (b.string[length(b.string)] != "/") {
      cdhit.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    cdhit.path <- ""
  }

  #Initial checks
  if (assembly.directory == output.directory){ stop("You should not overwrite the original contigs.") }
  if (is.null(target.file) == TRUE){ stop("A fasta file of targets to match to assembly contigs is needed.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  #Gets contig file names
  file.names = list.files(assembly.directory)

  #headers for the blast db
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

 # Sets up multiprocessing
  cl <- snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl <- floor(memory / threads)

  #Loop for cd-hit est reductions
  foreach::foreach(i = seq_along(file.names), .packages = c("foreach", "PhyloProcessR", "Biostrings", "data.table")) %dopar% {

  #Matching and processing for each sample
  #for (i in seq_along(file.names)) {

    #Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])
    species.dir = paste0(output.directory, "/", sample)

    #Creates species directory if none exists
    if (file.exists(species.dir) == F){ dir.create(species.dir) }

    #Checks if this has been done already
    if (overwrite == FALSE){
      if (file.exists(paste0(species.dir, "/", sample, "_matching-contigs.fa")) == T){
        print(paste0(sample, " already finished, skipping. Set overwrite to T if you want to overwrite."))
        next
      }
    }#end

    #########################################################################
    # Part A: reduce redundancy
    #########################################################################

    system(paste0(
      cdhit.path, "cd-hit-est -i ", assembly.directory, "/", file.names[i],
      " -o ", species.dir, "/", sample, "_red.fa -p 0 -T 1",
      " -n 5 -c 0.9 -M ", mem.cl * 100
    ))

    ### Read in data
    all.data = Biostrings::readDNAStringSet(file = paste0(species.dir, "/", sample, "_red.fa"), format = "fasta")

    names(all.data) = paste0("contig_", seq(seq_along(all.data)))

    # Writes the final loci
    final.loci = as.list(as.character(all.data))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(species.dir, "/", sample, "_rename.fa"),
      nbchar = 1000000, as.string = TRUE, open = "w"
    )

    #########################################################################
    #Part B: Blasting
    #########################################################################

    # Make blast database for the probe loci
    system(paste0(
      blast.path, "makeblastdb -in ", species.dir, "/", sample, "_rename.fa",
      " -parse_seqids -dbtype nucl -out ", species.dir, "/", sample, "_nucl-blast_db"
    ), ignore.stdout = quiet)

    # Matches samples to loci
    system(paste0(
      blast.path, "blastn -task dc-megablast -db ", species.dir, "/", sample, "_nucl-blast_db -evalue 0.001",
      " -query ", target.file, " -out ", species.dir, "/", sample, "_target-blast-match.txt",
      " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
      " -num_threads 1"
    ))

    #Need to load in transcriptome for each species and take the matching transcripts to the database
    system(paste0("rm ", species.dir, "/*nucl-blast_db*"))

    #Loads in match data
    match.data = data.table::fread(paste0(species.dir, "/", sample, "_target-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    data.table::setnames(match.data, headers)

    #Matches need to be greater than 12
    filt.data = match.data[match.data$matches > min.match.length,]
    #Percent identitiy must match 50% or greater
    filt.data = filt.data[filt.data$pident >= min.match.percent,]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      next }

    #Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    #Make sure the hit is greater than 50% of the reference length
    filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$qLen),]

    filt.data[filt.data$qName == "FrogCap_M24103",]

    #Reads in contigs
    contigs = all.data

    #########################################################################
    #Part B: Multiple sample contigs (tName) matching to one target (qName)
    #########################################################################
    #Pulls out
    contig.names = unique(filt.data[duplicated(filt.data$qName) == T,]$qName)

    #Saves non duplicated data
    good.data = filt.data[!filt.data$qName %in% contig.names,]

    #Only runs if there are duplicates
    fix.seq = Biostrings::DNAStringSet()
    if (length(contig.names) != 0){
      new.data = c()
      for (j in 1:length(contig.names)) {
        #Subsets data
        sub.match = filt.data[filt.data$qName %in% contig.names[j],]

        ########
        #Saves if they are on the same contig and same locus and fragmented for some reason
        ####################
        if (length(unique(sub.match$qName)) == 1 && length(unique(sub.match$tName)) == 1){
          new.qstart = min(sub.match$qStart, sub.match$qEnd)[1]
          new.qend = max(sub.match$qStart, sub.match$qEnd)[1]
          new.tstart = min(sub.match$tStart, sub.match$tEnd)[1]
          new.tend = max(sub.match$tStart, sub.match$tEnd)[1]
          sub.match$qStart = new.qstart
          sub.match$qEnd = new.qend
          sub.match$tStart = new.tstart
          sub.match$tEnd = new.tend
          sub.match$bitscore = sum(sub.match$bitscore)
          sub.match$matches = sum(sub.match$matches)
          new.data = rbind(new.data, sub.match[1,])
          next
        } #end if

        ########
        #Saves if they are two separate contigs but non-overlapping on the same locus; N repair
        ####################
        #Keep if they match to same contig, then not a paralog
        if (length(unique(sub.match$qName)) == 1){

          #Finds out if they are overlapping
          for (k in 1:nrow(sub.match)){
            new.start = min(sub.match$tStart[k], sub.match$tEnd[k])
            new.end = max(sub.match$tStart[k], sub.match$tEnd[k])
            sub.match$tStart[k] = new.start
            sub.match$tEnd[k] = new.end
          }#end k loop

          #If the number is negative then problem!
          hit.para = 0
          for (k in 1:(nrow(sub.match)-1)){
            if (sub.match$qStart[k+1]-sub.match$qEnd[k] < -30){ hit.para = 1 }
          }

          #If there are overlaps
          if (hit.para == 1){
            save.match = sub.match[sub.match$bitscore == max(sub.match$bitscore),]
            new.data = rbind(new.data, save.match)
            next
          }#end if

          #Adjacent and barely overlapping
          if (hit.para == 0){
            #Cuts the node apart and saves separately
            sub.match$qStart[1] = as.numeric(1)
            sub.match$tStart[1] = as.numeric(1)
            sub.match$qEnd[nrow(sub.match)] = sub.match$qLen[nrow(sub.match)]
            sub.match$tEnd[nrow(sub.match)] = sub.match$tLen[nrow(sub.match)]

            #Collects new sequence fragments
            spp.seq = contigs[names(contigs) %in% sub.match$tName]
            spp.seq = spp.seq[pmatch(sub.match$tName, names(spp.seq))]

            new.seq = Biostrings::DNAStringSet()
            for (k in 1:length(spp.seq)){
              n.pad = sub.match$qStart[k+1]-sub.match$qEnd[k]
              new.seq = append(new.seq, Biostrings::subseq(x = spp.seq[k], start = sub.match$tStart[k], end = sub.match$tEnd[k]) )
              if (is.na(n.pad) != T){ if (n.pad > 1){ new.seq = append(new.seq, Biostrings::DNAStringSet(paste0(rep("N", n.pad), collapse = "")) ) } }
            }#end kloop

            #Combine new sequence
            save.contig = Biostrings::DNAStringSet(paste0(as.character(new.seq), collapse = "") )
            names(save.contig) = paste0(sub.match$qName[1], "_:_", sub.match$tName[1], "_|_", sample)
            fix.seq = append(fix.seq, save.contig)
            next
          }#end if

          #Sets up the new contig location
          #Cuts the node apart and saves separately
          sub.match$tEnd = sub.match$tEnd+(sub.match$qSize-sub.match$qEnd)
          sub.contigs = contigs[names(contigs) %in% sub.match$qName]

          join.contigs<-DNAStringSet()
          for (k in 1:(nrow(sub.match)-1)){
            join.contigs<-append(join.contigs, sub.contigs[k])
            n.pad<-sub.match$tStart[k+1]-sub.match$tEnd[k]
            join.contigs<-append(join.contigs, DNAStringSet(paste(rep("N", n.pad), collapse = "", sep = "")) )
          }
          join.contigs<-append(join.contigs, sub.contigs[length(sub.contigs)])
          save.contig<-DNAStringSet(paste(as.character(join.contigs), collapse = "", sep = "") )

          #Saves final sequence
          names(save.contig)<-paste(loci.names[j], "_|_", sample, sep = "")
          fix.seq<-append(fix.seq, save.contig)


        }#end this if

        #Saves highest bitscore
        save.match = sub.match[sub.match$bitscore == max(sub.match$bitscore),]
        #Saves longest if equal bitscores
        save.match = save.match[abs(save.match$qStart-save.match$qEnd) == max(abs(save.match$qStart-save.match$qEnd)),]
        #saves top match here
        if (nrow(save.match) >= 2){  save.match = save.match[1,] }
        #Saves data
        new.data = rbind(new.data, save.match)
      } #end j

      #Saves final dataset
      save.data = rbind(good.data, new.data)
    } else { save.data = good.data }

    #########################################################################
    #Part D: Multiple targets (qName) matching to one sample contig (tName)
    #########################################################################

    #red.contigs = contigs[names(contigs) %in% filt.data$tName]
    dup.contigs = filt.data$tName[duplicated(filt.data$tName)]
    dup.match = filt.data[filt.data$tName %in% dup.contigs, ]
    dup.data = dup.match[order(dup.match$tName)]

    #Loops through each potential duplicate
    dup.loci = unique(dup.data$tName)

    if (length(dup.loci) != 0){
      for (j in 1:length(dup.loci)){
        #pulls out data that matches to multiple contigs
        sub.data = dup.data[dup.data$tName %in% dup.loci[j],]
        sub.data = sub.data[order(sub.data$tStart)]

        #Fixes direction and adds into data
        #Finds out if they are overlapping
        for (k in 1:nrow(sub.data)){
          new.start = min(sub.data$tStart[k], sub.data$tEnd[k])
          new.end = max(sub.data$tStart[k], sub.data$tEnd[k])
          sub.data$tStart[k] = new.start
          sub.data$tEnd[k] = new.end
        }#end k loop

        #Saves them if it is split up across the same locus
        if (length(unique(sub.data$tName)) == 1 && length(unique(sub.data$qName)) == 1){
          spp.seq = contigs[names(contigs) %in% sub.data$tName]
          names(spp.seq) = paste0(sub.data$qName[1], "_|_", sample)
          fix.seq = append(fix.seq, spp.seq)
          next
        }

        #Cuts the node apart and saves separately
        sub.data$tStart = sub.data$tStart-(sub.data$qStart-1)
        #If it ends up with a negative start
        sub.data$tStart[sub.data$tStart <= 0] = 1
        #Fixes ends
        sub.data$tEnd = sub.data$tEnd+(sub.data$qLen-sub.data$qEnd)

        #Fixes if the contig is smaller than the full target locus
        sub.data$tEnd[sub.data$tEnd >= sub.data$tLen] = sub.data$tLen[1]

        starts = c()
        ends = c()
        starts[1] = 1
        for (k in 1:(nrow(sub.data)-1)){
          ends[k] = sub.data$tEnd[k]+floor((sub.data$tStart[k+1]-sub.data$tEnd[k])/2)
          starts[k+1] = ends[k]+1
        } #end k loop
        ends = append(ends, sub.data$tLen[1])

        #Looks for overlapping contigs
        tmp = ends-starts
        if(length(tmp[tmp < 0 ]) != 0){
          sub.data = sub.data[sub.data$bitscore == max(sub.data$bitscore),]
          ends = sub.data$tEnd
          starts = sub.data$tStart
          # if (nrow(sub.data) != 1) { stop("ernor")}
        }

        #Collects new sequence fragments
        spp.seq = contigs[names(contigs) %in% sub.data$tName]
        new.seq = Biostrings::DNAStringSet()
        for (k in 1:length(starts)){ new.seq = append(new.seq, Biostrings::subseq(x = spp.seq, start = starts[k], end = ends[k]) ) }

        # #Sets up the new contig location
        # #Cuts the node apart and saves separately
        # sub.match$tEnd<-sub.match$tEnd+(sub.match$qSize-sub.match$qEnd)
        # sub.contigs<-contigs[names(contigs) %in% sub.match$qName]
        #
        # join.contigs<-DNAStringSet()
        # for (k in 1:(nrow(sub.match)-1)){
        #   join.contigs<-append(join.contigs, sub.contigs[k])
        #   n.pad<-sub.match$tStart[k+1]-sub.match$tEnd[k]
        #   join.contigs<-append(join.contigs, DNAStringSet(paste(rep("N", n.pad), collapse = "", sep = "")) )
        # }
        # join.contigs<-append(join.contigs, sub.contigs[length(sub.contigs)])
        # save.contig<-DNAStringSet(paste(as.character(join.contigs), collapse = "", sep = "") )

        #renames and saves
        names(new.seq) = paste0(sub.data$qName,"_|_", sample)
        fix.seq = append(fix.seq, new.seq)
      } #end j loop
    }#end if

    #Writes the base loci
    base.data = save.data[!save.data$qName %in% gsub("_\\:_.*", "", names(fix.seq)),]
    base.loci = contigs[names(contigs) %in% base.data$tName]
    sort.data = base.data[match(names(base.loci), base.data$tName),]
    #Name and finalize
    names(base.loci) = paste0(sort.data$qName, "_|_", sample)
    fin.loci = append(base.loci, fix.seq)
    fin.loci = fin.loci[Biostrings::width(fin.loci) >= min.match.length]

    #DUPES and numbers don't match up between contigs and table (dupes or not removed?)
    temp = fin.loci[duplicated(names(fin.loci)) == T]
    if(length(temp) != 0){
      dup.loci = fin.loci[names(fin.loci) %in% names(temp)]
      dup.loci = dup.loci[Biostrings::width(dup.loci) == max(Biostrings::width(dup.loci))][1]
      temp.fin = fin.loci[!names(fin.loci) %in% names(temp)]
      fin.loci = append(temp.fin, dup.loci)
    }#end duplicate if

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(fin.loci))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", sample, ".fa"), nbchar = 1000000, as.string = T
    )

    system(paste0("rm -r ", species.dir))

    print(paste0(sample, " target matching complete. ", length(final.loci), " targets found!"))

  }# end i loop

  snow::stopCluster(cl)

  ########################################################################
  # Output a single file for alignment
  ########################################################################

  #gets lists of directories and files with sample names
  file.names = list.files(assembly.directory)
  samples = gsub(".fa$", "", file.names)

  header.data = c("Sample", "startContigs", "annotatedContigs", "minLen", "maxLen", "meanLen")
  save.data = data.table::data.table(matrix(as.double(0), nrow = length(samples), ncol = length(header.data)))
  data.table::setnames(save.data, header.data)
  save.data[, Sample:=as.character(samples)]

  #Cycles through each assembly run and assesses each
  save.contigs = Biostrings::DNAStringSet()
  for (i in 1:length(samples)){

    #Sets up working directories for each species
    species.dir = paste0(output.directory, "/", samples[i])

    #Gets length of contigs
    og.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", samples[i], ".fa"))
    cd.contigs = Biostrings::readDNAStringSet(paste0(output.directory, "/", samples[i], ".fa"))

    #Gets the saved matching targets
    data.table::set(save.data, i =  match(samples[i], samples), j = match("Sample", header.data), value = samples[i] )
    data.table::set(save.data, i = match(samples[i], samples), j = match("startContigs", header.data), value = length(og.contigs) )
    data.table::set(save.data, i = match(samples[i], samples), j = match("annotatedContigs", header.data), value = length(cd.contigs) )

    data.table::set(save.data, i =  match(samples[i], samples), j = match("minLen", header.data), value = min(Biostrings::width(cd.contigs)) )
    data.table::set(save.data, i =  match(samples[i], samples), j = match("maxLen", header.data), value = max(Biostrings::width(cd.contigs)) )
    data.table::set(save.data, i =  match(samples[i], samples), j = match("meanLen", header.data), value = mean(Biostrings::width(cd.contigs)) )

    names(cd.contigs) = paste0(gsub("_:_.*", "", names(cd.contigs)), "_|_", samples[i])
    save.contigs = append(save.contigs, cd.contigs)

  }#End loop for things

  #Finds probes that match to two or more contigs
  final.loci = as.list(as.character(save.contigs))
  writeFasta(sequences = final.loci, names = names(final.loci),
             paste0(alignment.contig.name, "_to-align.fa"), nbchar = 1000000, as.string = T)

  #Saves combined, final dataset
  write.csv(save.data, file = paste0(alignment.contig.name, "_sample-assessment.csv"), row.names = F)

} #End function


#### ** IDEA
#### Make table of contig matches to tagets, coordinates in each and length

#### END SCRIPT
