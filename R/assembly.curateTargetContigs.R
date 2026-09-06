#' @title curateTargetContigs
#'
#' @description Curates the contigs of each sample against the target markers.
#'   The function reduces redundancy with cd-hit-est, matches the targets to the
#'   contigs with LAST, joins the fragments of one target that sit on separate
#'   contigs, and cuts apart a contig that spans more than one target. It writes
#'   one sequence per target for each sample, named after the target.
#'
#' @details This is the structural half of \code{annotateTargets}. It runs at the
#'   end of workflow 2, before variant calling, because the variant caller maps
#'   the reads back to these contigs. A contig that spans two targets collects
#'   the reads of both loci and gives wrong genotypes, and no later step can
#'   repair that. \code{annotateTargets} keeps the other half of the work: the
#'   paralog policy, the sample naming, and the alignment output.
#'
#'   Run this step with \code{retain.paralogs} in mind. This function does not
#'   drop a paralog copy, because the reads of a dropped copy map to the copy
#'   that remains and make false heterozygosity.
#'
#' @param assembly.directory path to a directory of per-sample contig files, one
#'   \code{<sample>.fa} file per sample.
#'
#' @param target.file path to the FASTA file of target markers.
#'
#' @param output.directory path to the directory for the curated contigs.
#'   Default: \code{"curated-contigs"}.
#'
#' @param min.match.percent minimum percent identity of a match. Default:
#'   \code{60}.
#'
#' @param min.match.length minimum alignment length in base pairs. Default:
#'   \code{60}.
#'
#' @param min.match.coverage minimum percentage of the target length that the
#'   matches must cover. The hits of one target are summed, and the test runs
#'   after the fragments are joined. A target split across two contigs therefore
#'   passes when the two fragments together cover enough of it. The default is
#'   permissive, because a later step can still remove a short locus, but this
#'   step cannot recover one it dropped. Default: \code{30}.
#'
#' @param search.method which program matches the target markers to the contigs.
#'   \code{"last"} (default) uses LAST, which matches a contig that is up to
#'   about 35 percent divergent from its target. \code{"blast"} uses
#'   \code{blastn dc-megablast}, which loses a contig past about 25 percent
#'   divergence. Numbers in HANDOFF.md.
#'
#' @param threads number of samples to process at the same time. Default:
#'   \code{1}.
#'
#' @param memory RAM in GB shared by the parallel samples. Default: \code{1}.
#'
#' @param blast.path path to the directory containing the BLAST executables.
#'   Only needed when \code{search.method = "blast"}. Default: \code{NULL}.
#'
#' @param last.path path to the directory containing \code{lastdb} and
#'   \code{lastal}. Only needed when \code{search.method = "last"}. Default:
#'   \code{NULL}.
#'
#' @param cdhit.path path to the directory containing \code{cd-hit-est}.
#'   Default: \code{NULL}.
#'
#' @param overwrite logical. \code{TRUE} runs a sample again when its output
#'   exists. Default: \code{FALSE}.
#'
#' @param quiet logical. \code{TRUE} hides the output of the external programs.
#'   Default: \code{TRUE}.
#'
#' @return Invisibly returns nothing. Writes one \code{<sample>.fa} file per
#'   sample to \code{output.directory}, with each sequence named after its
#'   target, and a match table per sample to
#'   \code{logs/sample_logs/<sample>_curation-matches.csv}.
#'
#' @export

curateTargetContigs = function(assembly.directory = NULL,
                               target.file = NULL,
                               output.directory = "curated-contigs",
                               min.match.percent = 60,
                               min.match.length = 60,
                               min.match.coverage = 30,
                               search.method = c("last", "blast"),
                               threads = 1,
                               memory = 1,
                               blast.path = NULL,
                               last.path = NULL,
                               cdhit.path = NULL,
                               overwrite = FALSE,
                               quiet = TRUE) {

  #Add the slash character to path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  if (is.null(cdhit.path) == FALSE){
    b.string = unlist(strsplit(cdhit.path, ""))
    if (b.string[length(b.string)] != "/") {
      cdhit.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { cdhit.path = "" }

  search.method = match.arg(search.method)
  last.path = .programPrefix(last.path)

  #Initial checks
  if (is.null(assembly.directory) == TRUE) { stop("Please provide a contig directory.") }
  if (is.null(target.file) == TRUE) { stop("Please provide a target file.") }
  if (file.exists(target.file) == FALSE) { stop("Target file not found.") }

  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists("logs/sample_logs")) {
    dir.create("logs/sample_logs", recursive = TRUE, showWarnings = FALSE)
  }

  #Gets contig file names
  file.names = list.files(assembly.directory)

  #headers for the search results
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
            "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  mem.cl <- floor(memory / threads)

  results = parallel::mclapply(seq_along(file.names), function(i) {
  tryCatch({

    #Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])

    #Checks if this has been done already (before creating any directory)
    if (overwrite == FALSE){
      if (file.exists(paste0(output.directory, "/", sample, ".fa")) == TRUE){
        print(paste0(sample, " already finished, skipping. Set overwrite = TRUE to redo."))
        return(NULL)
      }
    }#end

    # Temporary working directory kept inside logs so output.directory stays clean
    species.dir = paste0("logs/sample_logs/", sample)
    if (!dir.exists(species.dir)){ dir.create(species.dir, recursive = TRUE, showWarnings = FALSE) }

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

    # The database holds the contigs of this sample and the query is the target
    # file, so a hit names the target first and the contig second.
    # One thread per search, because the samples already run in parallel.
    if (search.method == "last") {
      .lastBuildDB(reference.file = paste0(species.dir, "/", sample, "_rename.fa"),
                   db.prefix = paste0(species.dir, "/", sample, "_last_db"),
                   lastdb.command = paste0(last.path, "lastdb"),
                   threads = 1,
                   quiet = quiet)

      .lastSearch(query.file = target.file,
                  db.prefix = paste0(species.dir, "/", sample, "_last_db"),
                  out.file = paste0(species.dir, "/", sample, "_target-blast-match.txt"),
                  lastal.command = paste0(last.path, "lastal"),
                  threads = 1,
                  quiet = quiet)
    } else {
      system(paste0(
        blast.path, "makeblastdb -in ", species.dir, "/", sample, "_rename.fa",
        " -parse_seqids -dbtype nucl -out ", species.dir, "/", sample, "_nucl-blast_db"
      ), ignore.stdout = quiet)

      system(paste0(
        blast.path, "blastn -task dc-megablast -db ", species.dir, "/", sample, "_nucl-blast_db -evalue 0.001",
        " -query ", target.file, " -out ", species.dir, "/", sample, "_target-blast-match.txt",
        " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
        " -num_threads 1"
      ))
    }

    # Remove the search database and the large intermediate contig files
    system(paste0("rm -f ", species.dir, "/*nucl-blast_db* ", species.dir, "/*_last_db*"))
    system(paste0("rm -f ",
      species.dir, "/", sample, "_red.fa ",
      species.dir, "/", sample, "_red.fa.clstr ",
      species.dir, "/", sample, "_rename.fa"
    ))

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

    # The coverage test is not applied here. It runs after Part C and Part D,
    # once the fragments of a target are joined. See the note there.

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
            spp.seq = spp.seq[match(sub.match$tName, names(spp.seq))]

            new.seq = Biostrings::DNAStringSet()
            for (k in 1:length(spp.seq)){
              n.pad = sub.match$qStart[k+1]-sub.match$qEnd[k]
              new.seq = append(new.seq, Biostrings::subseq(x = spp.seq[k], start = sub.match$tStart[k], end = sub.match$tEnd[k]) )
              if (is.na(n.pad) != T){ if (n.pad > 1){ new.seq = append(new.seq, Biostrings::DNAStringSet(paste0(rep("N", n.pad), collapse = "")) ) } }
            }#end kloop

            #Combine new sequence
            save.contig = Biostrings::DNAStringSet(paste0(as.character(new.seq), collapse = "") )
            names(save.contig) = sub.match$qName[1]
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
          names(spp.seq) = sub.data$qName[1]
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
        names(new.seq) = sub.data$qName
        fix.seq = append(fix.seq, new.seq)
      } #end j loop
    }#end if


    #########################################################################
    #Part E: Write the curated set
    #########################################################################

    # One sequence per target. The repaired sequences of Part C and the split
    # sequences of Part D replace the plain contig for those targets.
    fix.seq.final = append(fix.seq, fix.seq.para)
    base.data = save.data[!save.data$qName %in% names(fix.seq.final),]
    base.loci = contigs[names(contigs) %in% base.data$tName]
    sort.data = base.data[match(names(base.loci), base.data$tName),]
    names(base.loci) = sort.data$qName

    fin.loci = append(base.loci, fix.seq.final)
    fin.loci = fin.loci[Biostrings::width(fin.loci) >= min.match.length]

    # Coverage is summed over every hit of a target, and the test runs here, not
    # in Part B. A target that is split across two contigs has two hits that are
    # each too short on their own. A test before Part C drops both of them, and
    # the joined sequence that Part C would have made is lost.
    target.cover = tapply(filt.data$matches, filt.data$qName, sum)
    target.len   = tapply(filt.data$qLen, filt.data$qName, max)
    keep.targets = names(target.cover)[target.cover >=
                                       ((min.match.coverage / 100) * target.len)]
    fin.loci = fin.loci[names(fin.loci) %in% keep.targets]

    if (length(fin.loci) == 0) {
      print(paste0(sample, " had no curated contigs. Skipping"))
      return(NULL)
    }

    # A paralog keeps its own copy. Dropping one here would send its reads to the
    # copy that remains and make false heterozygosity in the variant caller.
    names(fin.loci) = make.unique(names(fin.loci), sep = "_")

    final.loci = as.list(as.character(fin.loci))
    PhyloProcessR::writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(output.directory, "/", sample, ".fa"), nbchar = 1000000, as.string = T
    )

    #------------------------------------------------------
    # Per-sample match log
    #------------------------------------------------------
    filt.log = as.data.frame(filt.data)[, c("qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
    filt.log$Sample = sample
    filt.log = filt.log[, c("Sample", "qName", "tName", "pident", "matches", "bitscore", "evalue", "qLen", "tLen")]
    write.csv(filt.log, file = paste0("logs/sample_logs/", sample, "_curation-matches.csv"), row.names = FALSE)

    print(paste0(sample, " curation complete. ", length(final.loci), " targets kept."))

    return(data.frame(
      Sample          = sample,
      DedupContigs    = length(all.data),
      TargetsMatched  = length(unique(filt.data$qName)),
      CuratedTargets  = length(final.loci),
      MeanPident      = round(mean(filt.data$pident), 2),
      MeanBitscore    = round(mean(filt.data$bitscore), 1),
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    print(paste0(file.names[i], " failed: ", conditionMessage(e)))
    return("failed")
  })
  }, mc.cores = threads)

  # A warning raised in a forked child never reaches the parent, so a failed
  # sample must be counted here.
  fail.count = sum(vapply(results, function(x) identical(x, "failed"), logical(1)))
  if (fail.count != 0) {
    print(paste0(fail.count, " of ", length(file.names),
                 " samples failed. See the messages above."))
  }

  #Writes the cross-sample summary
  summary.df = do.call(rbind, results[vapply(results, is.data.frame, logical(1))])
  if (is.null(summary.df) == FALSE && nrow(summary.df) != 0) {
    write.csv(summary.df, file = "logs/curateTargetContigs_summary.csv", row.names = FALSE)
  }

  return(invisible(NULL))
}#end function
