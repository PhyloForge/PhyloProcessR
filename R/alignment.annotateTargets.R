#' @title annotateTargets
#'
#' @description Annotates assembly contigs. It matches them to a set of target
#' marker sequences. For each sample, the function first removes redundant contigs
#' with CD-HIT-EST. It then searches the contigs against the target file.
#' \code{search.method} selects the program. The default is LAST. The function
#' filters the hits by percent identity, match length, and coverage. Contigs
#' that span multiple targets or targets that span multiple contigs are handled
#' by trimming or joining with N padding. The annotated contigs for each sample
#' are saved as a per-sample FASTA file in \code{output.directory}. A combined
#' FASTA file suitable for downstream alignment (named
#' \code{alignment.contig.name_to-align.fa}) and a summary CSV are written to
#' the working directory.
#'
#' @details The structural curation moved to \code{curateTargetContigs}, which
#'   runs at the end of workflow 2. That function joins the fragments of one
#'   target and cuts apart a contig that spans more than one target, before the
#'   variant caller maps reads to the contigs. This function keeps the paralog
#'   policy, the sample naming, and the alignment output. Give it contigs that
#'   \code{curateTargetContigs} has already curated.
#'
#' @param assembly.directory path to the directory containing per-sample contig FASTA files
#' (one file per sample, named \code{sampleName.fa}).
#'
#' @param target.file path to the FASTA file of target marker sequences used for the search
#' matching.
#'
#' @param alignment.contig.name base name (without extension) used for the combined output
#' FASTA and summary CSV files. Default "annotated-contigs-all".
#'
#' @param output.directory path to the directory where per-sample annotated contig files
#' will be saved. Default "annotated-contigs".
#'
#' @param min.match.percent minimum percent identity required to retain a hit.
#' Default 60.
#'
#' @param min.match.length minimum alignment length (in bp) required to retain a hit.
#' Default 60.
#'
#' @param min.match.coverage minimum proportion of the target sequence length that must be
#' covered by the hit (expressed as a percentage). Default 50.
#'
#' @param retain.paralogs logical. If TRUE, potential paralogs (multiple contigs matching
#' the same target) are retained by keeping the highest-bitscore hit. If FALSE, the
#' best-scoring contig is selected. Default FALSE.
#'
#' @param threads number of parallel threads to use. Default 1.
#'
#' @param memory total memory (in GB) to allocate across all threads. Default 1.
#'
#' @param search.method which program matches the target markers to the contigs.
#'   \code{"last"} (default) uses LAST, which matches a contig that is up to
#'   about 35 percent divergent from its target. \code{"blast"} uses the
#'   previous \code{blastn dc-megablast} search, which loses a contig past about
#'   25 percent divergence. Numbers in HANDOFF.md.
#'
#' @param last.path path to the directory containing \code{lastdb} and
#'   \code{lastal}. Only needed when \code{search.method = "last"}. If NULL the
#'   programs must be on the system PATH.
#'
#' @param blast.path path to the directory containing BLAST executables. Only
#'   needed when \code{search.method = "blast"}. If NULL, BLAST
#' tools are expected to be on the system PATH.
#'
#' @param cdhit.path path to the directory containing the CD-HIT-EST executable. If NULL,
#' CD-HIT-EST is expected to be on the system PATH.
#'
#' @param overwrite logical. If TRUE, previously completed samples are reprocessed;
#' if FALSE, they are skipped. Default FALSE.
#'
#' @param quiet logical. If TRUE, suppresses the search program screen output. Default TRUE.
#'
#' @return Writes per-sample annotated FASTA files to \code{output.directory}, a combined
#' FASTA file for alignment, and a summary CSV to the working directory. Two log files are
#' also written:
#' \itemize{
#'   \item \code{logs/sample_logs/<Sample>_blast-matches.csv} -- the filtered hit table
#'     for each sample (one row per hit: target, contig, pident, bitscore, evalue, lengths).
#'   \item \code{logs/annotateTargets_summary.csv} -- one row per sample summarising
#'     deduplicated contig count, number of targets matched, annotated target count, and
#'     mean/max identity and bitscore.
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
                            last.path = NULL,
                            search.method = c("last", "blast"),
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

  search.method = match.arg(search.method)
  last.path = .programPrefix(last.path)

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

    #Make sure the hit is greater than 50% of the reference length
    filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$qLen),]

    #Reads in contigs
    contigs = all.data

    #########################################################################
    #Part C: One contig per target
    #########################################################################
    # curateTargetContigs does the structural work at the end of workflow 2. It
    # joins the fragments of one target and cuts apart a contig that spans more
    # than one target, before the variant caller maps reads to the contigs. Only
    # the best match per target is needed here.
    # Bitscore first, then the longer match on a tie.
    data.table::setorderv(filt.data, c("qName", "bitscore", "matches"),
                          order = c(1L, -1L, -1L))
    save.data = filt.data[duplicated(filt.data$qName) == FALSE, ]

    # Part C and Part D of the old function produced these two sets. They stay,
    # empty, so the code below reads the same as it did.
    fix.seq = Biostrings::DNAStringSet()
    fix.seq.para = Biostrings::DNAStringSet()

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
    print(paste0(file.names[i], " failed: ", conditionMessage(e)))
    return("failed")
  })
  }, mc.cores = threads) # end i loop

  # A warning raised in a forked child never reaches the parent, so a failed
  # sample must be counted here.
  fail.count = sum(vapply(results, function(x) identical(x, "failed"), logical(1)))
  if (fail.count != 0) {
    print(paste0(fail.count, " of ", length(file.names),
                 " samples failed. See the messages above."))
  }

  ########################################################################
  # Write cross-sample summary log
  ########################################################################

  summary.df = do.call(rbind, results[vapply(results, is.data.frame, logical(1))])
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

    #Gets length of contigs
    og.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", samples[i], ".fa"))
    out.fa = paste0(output.directory, "/", samples[i], ".fa")
    if (!file.exists(out.fa)) {
      warning(samples[i], ": no annotated output file found -- sample had no targets passing filters, skipping.")
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
  write.csv(save.data, file = "logs/annotation_sample_summary.csv", row.names = F)

} #End function


#### ** IDEA
#### Make table of contig matches to tagets, coordinates in each and length

#### END SCRIPT
