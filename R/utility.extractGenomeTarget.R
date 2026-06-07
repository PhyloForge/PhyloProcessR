#' @title extractGenomeTarget
#'
#' @description Extracts genomic sequences from one or more genome assemblies
#'   that correspond to a set of target loci. Targets can be provided either as
#'   a FASTA file (matched via dc-megablast) or as a BED coordinate file.
#'   Matching genomic regions are extracted with optional flanking sequence and
#'   saved as per-genome FASTA files. Optionally also writes BED and match-table
#'   CSV outputs.
#'
#' @param genome.path path to a single genome FASTA file or to a directory of
#'   genome files to process.
#'
#' @param input.file path to the target sequence file; either a multi-sequence
#'   FASTA (when input.type = "fasta") or a BED coordinate file (when
#'   input.type = "bed").
#'
#' @param input.type character; "fasta" to BLAST targets against each genome,
#'   or "bed" to extract coordinates directly from a BED file.
#'
#' @param output.name name of the output directory where results are saved.
#'
#' @param output.bed logical; if TRUE a BED file of matched coordinates is
#'   written for each genome.
#'
#' @param bed.headers logical; if TRUE the BED file is treated as having a
#'   header row (only relevant when input.type = "bed").
#'
#' @param name.bed.names logical; reserved for naming BED entries from the
#'   target name column.
#'
#' @param output.table logical; if TRUE a CSV table of BLAST match statistics
#'   is written for each genome.
#'
#' @param match.by.chr logical; reserved for filtering matches to chromosome-
#'   level sequences only.
#'
#' @param duplicate.matches character; how to handle cases where multiple
#'   targets match the same genomic contig: "best" keeps the highest-scoring
#'   match, "all" keeps all matches, "none" discards duplicates.
#'
#' @param merge.matches logical; if TRUE, multiple BLAST hits from different
#'   parts of the same target to the same contig are merged into a single
#'   coordinate span before extraction.
#'
#' @param minimum.match.length integer minimum number of aligned bases required
#'   to retain a BLAST match.
#'
#' @param minimum.match.identity numeric minimum percent identity (0-100) for a
#'   BLAST match to be retained.
#'
#' @param minimum.match.coverage numeric minimum proportion of the target length
#'   that must be covered by the BLAST match.
#'
#' @param add.flanks integer number of base pairs to add upstream and downstream
#'   of each matched region when extracting sequences.
#'
#' @param genome.search.string optional character string to filter files in a
#'   genome directory by name (e.g. "_genomic.fna.gz").
#'
#' @param threads number of CPU threads to pass to blastn.
#'
#' @param memory amount of RAM in GB (currently reserved).
#'
#' @param overwrite logical; if TRUE existing output directories are deleted and
#'   recreated.
#'
#' @param quiet logical; if TRUE BLAST stdout is suppressed.
#'
#' @param blast.path system path to the directory containing blastn and
#'   makeblastdb; NULL searches the system PATH.
#'
#' @return invisibly; extracted FASTA sequences are written to per-genome
#'   sub-directories inside output.name.
#'
#' @export

extractGenomeTarget = function(genome.path = NULL,
                               input.file = NULL,
                               input.type = c("fasta", "bed"),
                               output.name = NULL,
                               output.bed = TRUE,
                               bed.headers = FALSE,
                               name.bed.names = TRUE,
                               output.table = TRUE,
                               match.by.chr = FALSE,
                               duplicate.matches = c("best", "all", "none"),
                               merge.matches = FALSE,
                               minimum.match.length = 100,
                               minimum.match.identity = 0.75,
                               minimum.match.coverage = 0.75,
                               add.flanks = 500,
                               genome.search.string = NULL,
                               threads = 1,
                               memory = 1,
                               overwrite = FALSE,
                               quiet = FALSE,
                               blast.path = NULL) {

  #todo: make single function for one sample
  #todo: blast type
  #todo: names

  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Microhylidae_SeqCap")
  # genome.path = "/Volumes/LaCie/Reference_Genomes/Amphibians/Frogs/Gastrophryne_carolinensis_GCA-027917425/GCA_027917425.1_aGasCar1.pri_genomic.fna.gz"
  # input.file = "COGEDA_81_marker.fasta"
  # #input.file = "/Volumes/Armored/Anolis_UCE/Anolis_genome/Anolis_carolinensis_final-table.bed"
  # output.name = "amphibians"
  # input.type = "fasta"
  # match.by.chr = TRUE
  # output.bed = TRUE
  # bed.headers = TRUE
  # output.table = TRUE
  # duplicate.matches = "best"
  # minimum.match.length = 50
  # minimum.match.identity = 0.50
  # minimum.match.coverage = 0.5
  # add.flanks = 500
  # threads = 8
  # memory = 24
  # overwrite = TRUE
  # quiet = FALSE
  # merge.matches = TRUE
  # blast.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # genome.search.string = "_genomic.fna.gz"
#
  # setwd("/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/2_sea_snake_genomes")
  # genome.path = "/Users/chutter/Dropbox/Research/1_Main-Projects/1_Collaborative-Projects/Elapid_Probes/Gene_based/genomes/GCA_019472885_1_HCur_v2_genomic.fa"
  # input.file = "Hcur_coordinates.txt"
  # blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin/"
  #
  # output.name = "Hcur"
  # input.type = "bed"
  # output.bed = TRUE
  # bed.headers = FALSE
  # output.table = FALSE
  # match.by.chr = FALSE
  # duplicate.matches = "none"
  # merge.matches = FALSE
  # minimum.match.length = 100
  # minimum.match.identity = 0.75
  # add.flanks = 30
  # genome.search.string = "_genomic.fna.gz"
  # threads = 8
  # memory = 16
  # overwrite = FALSE
  # quiet = FALSE

  #Same adds to blast path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  # Validate choice parameters
  input.type        = match.arg(input.type)
  duplicate.matches = match.arg(duplicate.matches)

  #Initial checks
  if (output.name == genome.path){ stop("You should not overwrite the original genome file.") }
  if (output.name == input.file){ stop("You should not overwrite the original target file.") }
  if (is.null(input.file) == T){ stop("A fasta file of targets to match to assembly contigs is needed.") }
  if (is.null(genome.path) == T){ stop("A fasta file of targets to match to assembly contigs is needed.") }

  if (dir.exists(output.name) == TRUE) {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.name))
      dir.create(output.name)
    }
  } else { dir.create(output.name) }

  #Gets file names if they are directory or a single file
  if (dir.exists(genome.path) == TRUE) {
    file.names = list.files(genome.path, recursive = T)
    genome.files = file.names[grep("\\.fai$", file.names, invert = T)]
    if (!is.null(genome.search.string)){
      genome.files = genome.files[grep(genome.search.string, genome.files)]
    }
  } else {
    genome.files = genome.path
    genome.path = dirname(genome.files)
    genome.files = basename(genome.files)
  }#end else


  # if (match.by.chr == TRUE){
  #   chrom.samples = file.names[grep("chromosomes", file.names)]
  #   chrom.names = gsub("/chromosomes.*", "",  chrom.samples)
  #   chrom.names = unique(gsub(".*/", "", chrom.names))
  #   #remove and add from original files
  #   genome.files = genome.files[grep(paste(chrom.names, collapse = "|"), genome.files, invert = T)]
  #   genome.files = append(genome.files, chrom.samples)
  # }#end if


  if (input.type == "fasta"){

    #headers for the blast db
    headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

    ### Reads in the additional stuff
    target.markers = Biostrings::readDNAStringSet(input.file)  # loads up fasta file

    #Make blast database for the probe loci
    system(paste0(blast.path, "makeblastdb -in ", input.file,
                  " -parse_seqids -dbtype nucl -out ", output.name, "/target_nucl-blast_db"), ignore.stdout = quiet)
  }#end fasta


  #Loops through each locus and does operations on them
  for (i in 1:length(genome.files)) {
    #Sets up working directories for each species
    sample = gsub(pattern = "\\.fna\\.gz$|\\.fa\\.gz$|\\.fna$|\\.fa$|\\.gz$", replacement = "", x = genome.files[i])
    sample = gsub(pattern = ".*/", replacement = "", x = sample)
    species.dir = paste0(output.name, "/", sample)

    #Creates species directory if none exists
    if (dir.exists(species.dir) == F){ dir.create(species.dir) }

    if (input.type == "fasta"){

      #Checks if this has been done already
      if (overwrite == FALSE){
        if (file.exists(paste0(species.dir, "/", sample, "_target-matches.fa")) == T){ next }
      }#end

      #########################################################################
      #Part A: Blasting
      #########################################################################

      zipped.up = grep(".gz$", genome.files[i])

      if (length(zipped.up) == 1){
        system(paste0("gzip -dc ", genome.path, "/", genome.files[i], " > ",
                                  species.dir, "/", sample, "_genome.fa") )
        species.genome.path = paste0(species.dir, "/", sample, "_genome.fa")
        } else {
        species.genome.path = paste0(genome.path, "/", genome.files[i])
      }#end else

      #Matches samples to loci
      system(paste0(blast.path, "blastn -task dc-megablast -db ", output.name, "/target_nucl-blast_db -evalue 0.001",
                    " -query ", species.genome.path,
                    " -out ", species.dir, "/", sample, "_target-blast-match.txt",
                    " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                    " -num_threads ", threads))

      #Matches samples to proteins
      # system(paste0("tblastn -task tblastn-fast -db ", sample, "_nucl-blast_db -evalue 0.001 -seg no",
      #               " -query ", prot.file, " -out ", sample, "_prot-match.txt",
      #               " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
      #               " -num_threads ", threads))
      #  #Matches need to be greater than 12
      #filt.data = match.data[match.data$matches > 12,]
      #Percent identitiy must match 50% or greater
      #filt.data = filt.data[filt.data$pident >= 50,]

      #Loads in match data
      match.data = data.table::fread(paste0(species.dir, "/", sample, "_target-blast-match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

      #Guard: empty BLAST output
      if (nrow(match.data) == 0) {
        print(paste0(sample, " had no matches. Skipping"))
        next
      }
      data.table::setnames(match.data, headers)

      #Matches need to be >= minimum.match.length bp
      filt.data = match.data[match.data$matches >= minimum.match.length,]
      #Remove identical rows
      filt.data = filt.data[duplicated(filt.data) != T,]
      #Percent identity: minimum.match.identity is 0-1 fraction; pident is 0-100
      filt.data = filt.data[filt.data$pident >= minimum.match.identity * 100,]

      #Make sure the hit is greater than 50% of the reference length
      filt.data = filt.data[filt.data$matches >= ( (minimum.match.coverage) * filt.data$tLen),]

      #Skips if no matches
      if (nrow(filt.data) == 0) {
        print(paste0(sample, " had no matches. Skipping"))
        next }

      #Sorting: exon name, contig name, bitscore higher first, evalue
      data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

      #########################################################################
      #Part B: Merge tName that match to multiple same contig qName
      #########################################################################

      if (merge.matches == TRUE){

        merge.names = filt.data$tName[duplicated(filt.data$tName)]
        merge.match = filt.data[filt.data$tName %in% merge.names,]
        merge.data = merge.match[order(merge.match$tName)]

        #Loops through each potential duplicate
        merge.loci = unique(merge.data$tName)

        if (length(merge.loci) != 0){
          save.merge = c()
          del.merge = c()
          for (j in 1:length(merge.loci)){
            #pulls out data that matches to multiple contigs
            sub.data = merge.data[merge.data$tName %in% merge.loci[j],]
            sub.data = sub.data[order(sub.data$tStart)]

            dup.merge = unique(sub.data$qName)
            del.merge = append(del.merge, merge.loci[j])

            for (k in 1:length(dup.merge)){
              #merge duplicates
              m.data = sub.data[sub.data$qName %in% dup.merge[k],]
              m.data = m.data[order(m.data$tStart)]

              if (nrow(m.data) == 1){ save.merge = rbind(save.merge, m.data) }

              if (nrow(m.data) != 1){
                #Merge into one record
                m.data$tStart[1] = min(c(m.data$tEnd, m.data$tStart))
                m.data$tEnd[1] = max(c(m.data$tEnd, m.data$tStart))
                m.data$qStart[1] = min(c(m.data$qEnd, m.data$qStart))
                m.data$qEnd[1] = max(c(m.data$qEnd, m.data$qStart))
                m.data$bitscore[1] = sum(m.data$bitscore)
                m.data$matches[1] = sum(m.data$matches)
                m.data$misMatches[1] = sum(m.data$misMatches)
                m.data$gapopen[1] = sum(m.data$gapopen)

                save.merge = rbind(save.merge, m.data[1,])
              }#end if

            }#end k loop

          }#end j

          #Filters to new dataset
          filt.data = filt.data[!filt.data$tName %in% del.merge,]
          filt.data = rbind(filt.data, save.merge)
        }#end if
      }#end merge.matches if

      #########################################################################
      #Part C: Multiple targets (qName) matching to one sample contig (tName)
      #########################################################################

      dup.names = filt.data$tName[duplicated(filt.data$tName)]
      dup.match = filt.data[filt.data$tName %in% dup.names,]
      dup.data = dup.match[order(dup.match$tName)]

      #Loops through each potential duplicate
      dup.loci = unique(dup.data$tName)

      if (length(dup.loci) != 0){
        save.dup = c()
        for (j in 1:length(dup.loci)){
          #pulls out data that matches to multiple contigs
          sub.data = dup.data[dup.data$tName %in% dup.loci[j],]
          sub.data = sub.data[order(sub.data$tStart)]

          #Here decide what to do with the duplicate matches based on selection
          if (duplicate.matches == "none"){ next }

          #Saves all the decent matches
          if (duplicate.matches == "all"){
            save.temp = sub.data
          }

          #Here find the best match for this set of duplicates
          if (duplicate.matches == "best"){
            save.temp = sub.data[sub.data$bitscore == max(sub.data$bitscore),][1,]
          }

          save.dup = rbind(save.dup, save.temp)
        } #end j loop

        temp.filt = filt.data[!filt.data$tName %in% dup.names,]
        if (duplicate.matches == "none") {
          filt.data = temp.filt
        } else {
          filt.data = rbind(temp.filt, save.dup)
        }

      }#end if

      final.table = filt.data[order(filt.data$tName)]

      #Extracts the genomic data using the final.table coordinates
      Rsamtools::indexFa(species.genome.path)
      fa = Rsamtools::FaFile(species.genome.path)

      #adds genome columns — initialise to NA so skipped rows are filtered below
      final.table[, gStart:=NA_real_]
      final.table[, gEnd:=NA_real_]
      header.data = colnames(final.table)

      for (j in 1:nrow(final.table)){
        #subsets data
        sub.match = final.table[j,]

        #Gets starts and stops
        gen.start = min(c(sub.match$qStart, sub.match$qEnd)) - add.flanks
        gen.end  = max(c(sub.match$qEnd, sub.match$qStart)) + add.flanks
        if (gen.start <= 0){ gen.start = as.numeric(1) }
        if (gen.end > sub.match$qLen) { gen.end = sub.match$qLen }
        # Skip if the extracted span is implausibly large (> 100x target length)
        # Leave gStart/gEnd as NA so this row is removed by the filter below
        if (gen.end - gen.start >= max(sub.match$tLen) * 100){ next }

        ## Save in table here
        data.table::set(final.table, i =  as.integer(j),
                        j = match("gStart", header.data), value = gen.start )

        data.table::set(final.table, i =  as.integer(j),
                        j = match("gEnd", header.data), value = gen.end )

      }#end j loop

      #removes rows that were skipped (gStart is still NA)
      final.table = final.table[!is.na(final.table$gStart),]

      #gets ranges of stuff and obtains sequences
      gr.ranges = data.frame(seqnames = final.table$qName,
                             start = final.table$gStart,
                             end = final.table$gEnd)
      fin.ranges = GenomicRanges::makeGRangesFromDataFrame(gr.ranges)

      #gets ranges of stuff and obtains sequences
      target.seqs = BSgenome::getSeq(fa, fin.ranges)
      names(target.seqs) = paste0(final.table$tName, "_|_", sample)

      #Writes the table
      if (output.table == TRUE){
        write.csv(final.table, file = paste0(species.dir, "/", sample, "_final-table.csv"),
                  row.names = FALSE)
      }

      #Writes the bed
      if (output.bed == TRUE){
        bam.table = data.frame(chr = final.table$qName, start = final.table$gStart, stop = final.table$gEnd)
        write.table(bam.table, file = paste0(species.dir, "/", sample, "_final-table.bed"),
                    row.names = FALSE,
                    quote = F,
                    sep = "\t")
      }#end if


    }#end input type fasta

    if (input.type == "bed"){

      zipped.up = grep(".gz$", genome.files[i])

      if (length(zipped.up) == 1){
        system(paste0("gzip -dc ", genome.path, "/", genome.files[i], " > ",
                      species.dir, "/", sample, "_genome.fa") )
        species.genome.path = paste0(species.dir, "/", sample, "_genome.fa")
      } else {
        species.genome.path = paste0(genome.path, "/", genome.files[i])
      }#end else

      #Extracts the genomic data using the final.table coordinates
      Rsamtools::indexFa(species.genome.path)
      fa = Rsamtools::FaFile(species.genome.path)

      #gets ranges of stuff and obtains sequences
      final.table = read.table(input.file, header = bed.headers)

      if (bed.headers == TRUE){
        gr.ranges = data.frame(seqnames = final.table$chrom,
                               start = as.integer(final.table$start),
                               end = as.integer(final.table$stop) )
        fin.ranges = GenomicRanges::makeGRangesFromDataFrame(gr.ranges)

        #gets ranges of stuff and obtains sequences
        target.seqs = BSgenome::getSeq(fa, fin.ranges)
        names(target.seqs) = paste0(final.table$name, "_|_", sample)

      } else {
        gr.ranges = data.frame(seqnames = final.table[,1],
                               start = as.integer(final.table[,2]),
                               end = as.integer(final.table[,3]) )
        fin.ranges = GenomicRanges::makeGRangesFromDataFrame(gr.ranges)

        #gets ranges of stuff and obtains sequences
        target.seqs = BSgenome::getSeq(fa, fin.ranges)
        names(target.seqs) = paste0(gr.ranges$seqnames, "_", gr.ranges$start, "-",gr.ranges$end, "_|_", sample)
      }#end else

    }#end bed

    #Writes the final loci
    final.loci = as.list(as.character(target.seqs))
    PhyloProcessR::writeFasta(sequences = final.loci, names = names(final.loci),
               paste0(species.dir, "/", sample, "_target-matches.fa"),
               nbchar = 1000000, as.string = T)

    print(paste0(sample, " Matching to the genome complete. ", length(final.loci), " targets extracted!"))

    if (length(zipped.up) == 1){
      system(paste0("rm ", species.dir, "/", sample, "_genome.fa"))
      fai.path = paste0(species.dir, "/", sample, "_genome.fa.fai")
      if (file.exists(fai.path)) { system(paste0("rm ", fai.path)) }
    }

  }# end i loop

} #End function


#END SCRIPT
