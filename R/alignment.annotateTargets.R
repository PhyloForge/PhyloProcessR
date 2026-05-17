#' @title annotateTargets
#'
#' @description Annotates assembly contigs by matching them to a set of target marker
#' sequences using BLAST. For each sample, contigs are first deduplicated with CD-HIT-EST,
#' then BLASTed against the target file. BLAST hits are filtered by percent identity,
#' match length, and coverage. Contigs that span multiple targets or targets that span
#' multiple contigs are handled by trimming or joining with N padding. The annotated
#' contigs for each sample are saved as a per-sample FASTA file in \code{output.directory}.
#' A combined FASTA file suitable for downstream alignment (named
#' \code{alignment.contig.name_to-align.fa}) and a summary CSV are written to the working
#' directory.
#'
#' @param assembly.directory path to the directory containing per-sample contig FASTA files
#' (one file per sample, named \code{sampleName.fa}).
#'
#' @param target.file path to the FASTA file of target marker sequences used for BLAST
#' matching.
#'
#' @param alignment.contig.name base name (without extension) used for the combined output
#' FASTA and summary CSV files. Default "annotated-contigs-all".
#'
#' @param output.directory path to the directory where per-sample annotated contig files
#' will be saved. Default "annotated-contigs".
#'
#' @param min.match.percent minimum BLAST percent identity required to retain a hit.
#' Default 60.
#'
#' @param min.match.length minimum BLAST alignment length (in bp) required to retain a hit.
#' Default 60.
#'
#' @param min.match.coverage minimum proportion of the target sequence length that must be
#' covered by the BLAST hit (expressed as a percentage). Default 50.
#'
#' @param retain.paralogs logical. If TRUE, potential paralogs (multiple contigs matching
#' the same target) are retained by keeping the highest-bitscore hit. If FALSE, the
#' best-scoring contig is selected. Default FALSE.
#'
#' @param threads number of parallel threads to use. Default 1.
#'
#' @param memory total memory (in GB) to allocate across all threads. Default 1.
#'
#' @param blast.path path to the directory containing BLAST executables. If NULL, BLAST
#' tools are expected to be on the system PATH.
#'
#' @param cdhit.path path to the directory containing the CD-HIT-EST executable. If NULL,
#' CD-HIT-EST is expected to be on the system PATH.
#'
#' @param overwrite logical. If TRUE, previously completed samples are reprocessed;
#' if FALSE, they are skipped. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses BLAST screen output. Default TRUE.
#'
#' @return Writes per-sample annotated FASTA files to \code{output.directory}, a combined
#' FASTA file for alignment, and a summary CSV to the working directory. Two log files are
#' also written:
#' \itemize{
#'   \item \code{logs/sample_logs/<Sample>_blast-matches.csv} — the filtered BLAST table
#'     for each sample (one row per hit: target, contig, pident, bitscore, evalue, lengths).
#'   \item \code{logs/annotateTargets_summary.csv} — one row per sample summarising
#'     deduplicated contig count, number of targets matched, annotated target count, and
#'     mean/max BLAST identity and bitscore.
#' }
#' No value is returned to R.
#'
#' @export

annotateTargets = function(assembly.directory = NULL,
                            target.file = NULL,
                            alignment.contig.name = "annotated-contigs-all",
                            output.directory = "annotated-contigs",
                            min.match.percent = 60,
                            min.match.length = 60,
                            min.match.coverage = 50,
                            retain.paralogs = FALSE,
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

#
#   blast.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
#   cdhit.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
#
#
#   # #Debug setup
#   setwd("/Users/chutter/Downloads")
#   assembly.directory <- "/Users/chutter/Downloads/5_iupac-contigs"
#   target.file = "/Users/chutter/Dropbox/Research/1_Main-Projects/0_Working-Projects/Anura_Phylogeny/Final_Files/FINAL_marker-seqs_May20-2023.fa"
#   output.directory = "annotated_paralogs"
#
#   #
#   # #Main settings
#   threads = 4
#   memory = 20
#   trim.target = FALSE
#   overwrite = FALSE
#   quiet = TRUE
#   retain.paralogs = FALSE
#
#   # #tweak settings (make some statements to check these)
#   min.match.percent = 60
#   min.match.length = 70
#   min.match.coverage = 35
  #

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
  if (file.exists(target.file) == FALSE){ stop("Target file not found. Please check path / use full path.") }

  if (dir.exists(output.directory) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } else { dir.create(output.directory) }

  if (!dir.exists("logs/sample_logs")) {
    dir.create("logs/sample_logs", recursive = TRUE, showWarnings = FALSE)
  }

  #Gets contig file names
  file.names = list.files(assembly.directory)

  #headers for the blast db
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  mem.cl <- floor(memory / threads)

  #Loop for cd-hit est reductions
  results = parallel::mclapply(seq_along(file.names), function(i) {
  tryCatch({

    #Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])
    species.dir = paste0(output.directory, "/", sample)

    #Creates species directory if none exists
    if (file.exists(species.dir) == F){ dir.create(species.dir) }

    #Checks if this has been done already
    if (overwrite == FALSE){
      if (file.exists(paste0(output.directory, "/", sample, ".fa")) == TRUE){
        print(paste0(sample, " already finished, skipping. Set overwrite = TRUE to redo."))
        return(NULL)
      }
    }#end

    #########################################################################
    # Part A: reduce redundancy
    #########################################################################

    system(paste0(
      cdhit.path, "cd-hit-est -i ", assembly.directory, "/", file.names[i],
      " -o ", species.dir, "/", sample, "_red.fa -p 0 -T 1",
      " -n 8 -c 0.9 -M ", mem.cl * 1000
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
      return(NULL)
      }

    #Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    #Make sure the hit is greater than 50% of the reference length
    filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$qLen),]

    #Reads in contigs
    contigs = all.data

    #########################################################################
    #Part C: Multiple sample contigs (tName) matching to one target (qName)
    #########################################################################
    #Pulls out
    target.names = unique(filt.data[duplicated(filt.data$qName) == T,]$qName)

    #Saves non duplicated data
    good.data = filt.data[!filt.data$qName %in% target.names,]

    #Only runs if there are duplicates
    fix.seq = Biostrings::DNAStringSet()
    if (length(target.names) != 0){
      new.data = c()
      for (j in 1:length(target.names)) {
        #Subsets data
        sub.match = filt.data[filt.data$qName %in% target.names[j],]

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

    fix.seq.para = fix.seq

    #########################################################################
    #Part D: Multiple targets (qName) matching to one sample contig (tName)
    #########################################################################

    #red.contigs = contigs[names(contigs) %in% filt.data$tName]
    dup.contigs = filt.data$tName[duplicated(filt.data$tName)]
    dup.match = filt.data[filt.data$tName %in% dup.contigs, ]
    dup.data = dup.match[order(dup.match$tName)]

    #Loops through each potential duplicate
    dup.loci = unique(dup.data$tName)

    fix.seq = Biostrings::DNAStringSet()
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

    #########################################################################
    #Part E: Keep paralogs or no
    #########################################################################

    #Keeps potential paralogs
    if (retain.paralogs == TRUE){

      target.names = unique(filt.data[duplicated(filt.data$qName) == T,]$qName)

      #Saves non duplicated data
      good.data = filt.data[!filt.data$qName %in% target.names,]

      save.data = c()
      for (j in 1:length(target.names)){

        temp.data = filt.data[filt.data$qName %in% target.names[j],]
        temp.save = temp.data[temp.data$bitscore == max(temp.data$bitscore),][1,]
        save.data = rbind(save.data, temp.save)
      }

      #Name and finalize
      comb.data = rbind(good.data, save.data)
      #base.loci = contigs[match(base.data$tName, names(contigs))]
      #names(base.loci) = paste0(base.data$qName, "_|_", sample)

      #fin.loci = append(base.loci, fix.seq)
      #fin.loci = fin.loci[Biostrings::width(fin.loci) >= min.match.length]
      #sort.data = base.data[match(names(base.loci), base.data$tName),]


      base.data = comb.data[!comb.data$qName %in% gsub("_\\|_.*", "", names(fix.seq)),]
      base.loci = contigs[names(contigs) %in% base.data$tName]
      sort.data = base.data[match(names(base.loci), base.data$tName),]
      #Name and finalize
      names(base.loci) = paste0(sort.data$qName, "_|_", sample)
      fin.loci = append(base.loci, fix.seq)
      fin.loci = fin.loci[Biostrings::width(fin.loci) >= min.match.length]

      #DUPES and numbers don't match up between contigs and table (dupes or not removed?)
      temp = fin.loci[duplicated(names(fin.loci)) == T]
      if(length(temp) != 0){

        dup.names = unique(names(temp))
        save.temp = Biostrings::DNAStringSet()
        for (j in 1:length(dup.names)){

          temp.data = fin.loci[names(fin.loci) %in% dup.names[j]]
          best.temp = temp.data[Biostrings::width(temp.data) == max(Biostrings::width(temp.data))][1]
          save.temp = append(save.temp, best.temp)
        }# end j loop

        temp.fin = fin.loci[!names(fin.loci) %in% names(temp)]
        fin.loci = append(temp.fin, save.temp)
      }#end duplicate if

      #Finds probes that match to two or more contigs
      final.loci = as.list(as.character(fin.loci))
      PhyloProcessR::writeFasta(
        sequences = final.loci, names = names(final.loci),
        paste0(output.directory, "/", sample, ".fa"), nbchar = 1000000, as.string = T
      )

      #------------------------------------------------------
      # Per-sample blast log
      #------------------------------------------------------
      filt.log = as.data.frame(filt.data)[, c("qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
      filt.log$Sample = sample
      filt.log = filt.log[, c("Sample", "qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
      write.csv(filt.log, file = paste0("logs/sample_logs/", sample, "_blast-matches.csv"), row.names = FALSE)

      system(paste0("rm -r ", species.dir))

      print(paste0(sample, " target matching complete. ", length(final.loci), " targets found!"))

      return(data.frame(
        Sample           = sample,
        DedupContigs     = length(all.data),
        TargetsMatched   = length(unique(filt.data$qName)),
        AnnotatedTargets = length(final.loci),
        MeanPident       = round(mean(filt.data$pident), 2),
        MeanBitscore     = round(mean(filt.data$bitscore), 1),
        MaxBitscore      = round(max(filt.data$bitscore), 1),
        stringsAsFactors = FALSE
      ))

    }#end if statement

    #Writes the base loci
    fix.seq.final = append(fix.seq, fix.seq.para)
    base.data = save.data[!save.data$qName %in% gsub("_\\|_.*", "", names(fix.seq.final)),]
    base.loci = contigs[names(contigs) %in% base.data$tName]
    sort.data = base.data[match(names(base.loci), base.data$tName),]
    #Name and finalize
    names(base.loci) = paste0(sort.data$qName, "_|_", sample)
    fin.loci = append(base.loci, fix.seq.final)
    fin.loci = fin.loci[Biostrings::width(fin.loci) >= min.match.length]

    #DUPES and numbers don't match up between contigs and table (dupes or not removed?)
    temp = fin.loci[duplicated(names(fin.loci)) == T]
    if(length(temp) != 0){
      dup.names = unique(names(temp))
      save.temp = Biostrings::DNAStringSet()
      for (j in 1:length(dup.names)){
        temp.data = fin.loci[names(fin.loci) %in% dup.names[j]]
        best.temp = temp.data[Biostrings::width(temp.data) == max(Biostrings::width(temp.data))][1]
        save.temp = append(save.temp, best.temp)
      }# end j loop
      temp.fin = fin.loci[!names(fin.loci) %in% names(temp)]
      fin.loci = append(temp.fin, save.temp)
    }#end duplicate if

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(fin.loci))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", sample, ".fa"), nbchar = 1000000, as.string = T
    )

    #------------------------------------------------------
    # Per-sample blast log
    #------------------------------------------------------
    filt.log = as.data.frame(filt.data)[, c("qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
    filt.log$Sample = sample
    filt.log = filt.log[, c("Sample", "qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
    write.csv(filt.log, file = paste0("logs/sample_logs/", sample, "_blast-matches.csv"), row.names = FALSE)

    system(paste0("rm -r ", species.dir))

    print(paste0(sample, " target matching complete. ", length(final.loci), " targets found!"))

    return(data.frame(
      Sample           = sample,
      DedupContigs     = length(all.data),
      TargetsMatched   = length(unique(filt.data$qName)),
      AnnotatedTargets = length(final.loci),
      MeanPident       = round(mean(filt.data$pident), 2),
      MeanBitscore     = round(mean(filt.data$bitscore), 1),
      MaxBitscore      = round(max(filt.data$bitscore), 1),
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    warning(file.names[i], " failed: ", conditionMessage(e))
  })
  }, mc.cores = threads) # end i loop

  ########################################################################
  # Write cross-sample summary log
  ########################################################################

  summary.df = do.call(rbind, results[!sapply(results, is.null)])
  if (!is.null(summary.df) && nrow(summary.df) > 0) {
    out.csv = "logs/annotateTargets_summary.csv"
    if (file.exists(out.csv)) {
      existing = read.csv(out.csv, stringsAsFactors = FALSE)
      existing = existing[!existing$Sample %in% summary.df$Sample, ]
      summary.df = rbind(existing, summary.df)
    }
    write.csv(summary.df, file = out.csv, row.names = FALSE)
  }

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
    out.fa = paste0(output.directory, "/", samples[i], ".fa")
    if (!file.exists(out.fa)) {
      warning(samples[i], ": no annotated output file found — sample had no targets passing filters, skipping.")
      next
    }
    cd.contigs = Biostrings::readDNAStringSet(out.fa)

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
  PhyloProcessR::writeFasta(sequences = final.loci, names = names(final.loci),
             paste0(alignment.contig.name, "_to-align.fa"), nbchar = 1000000, as.string = T)

  #Saves combined, final dataset
  write.csv(save.data, file = paste0(alignment.contig.name, "_sample-assessment.csv"), row.names = F)

} #End function


#### ** IDEA
#### Make table of contig matches to tagets, coordinates in each and length

#### END SCRIPT
