#' @title extractGenomeTarget
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

  # setwd("/Volumes/LaCie/Ultimate_FrogCap")
  # genome.path = "/Volumes/LaCie/Reference_Genomes"
  # input.file = "/Volumes/LaCie/Anolis/A_carolinensus_RELEC.fa"
  # #input.file = "/Volumes/Armored/Anolis_UCE/Anolis_genome/Anolis_carolinensis_final-table.bed"
  # output.name = "amphibians"
  # input.type = "fasta"
  # match.by.chr = TRUE
  # output.bed = TRUE
  # bed.headers = TRUE
  # output.table = TRUE
  # duplicate.matches = "best"
  # minimum.match.length = 100
  # minimum.match.identity = 0.60
  # add.flanks = 500
  # threads = 8
  # memory = 24
  # overwrite = FALSE
  # quiet = FALSE
  # merge.matches = TRUE
  # blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # genome.search.string = "_genomic.fna.gz"

  # setwd("/Users/chutter/Dropbox/Elapid_probe_design/Gene_based")
  # genome.path = "/Users/chutter/Dropbox/Elapid_probe_design/Gene_based/genomes/GCA_019473425_1_HCya_v2_genomic.fa"
  # input.file = "Hcya_coordinates.txt"
  # blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin/"
  # output.name = "Hcur"
  # input.type = "bed"
  # match.by.chr = TRUE
  # output.bed = TRUE
  # bed.headers = TRUE
  # output.table = TRUE
  # duplicate.matches = "best"
  # minimum.match.length = 100
  # minimum.match.identity = 0.60
  # add.flanks = 500
  # threads = 8
  # memory = 24
  # overwrite = FALSE
  # quiet = FALSE
  # merge.matches = TRUE
  # blast.path = "/Users/chutter/Bioinformatics/conda-envs/PhyloCap/bin"
  # genome.search.string = "_genomic.fna.gz"
#
#   genome.path = genome.dir
#   input.file = paste0(work.dir, "/BUSCO_all-reference.fa")
#   input.type = "fasta"
#   output.name = out.dir
#   output.bed = TRUE
#   bed.headers = FALSE
#   name.bed.names = TRUE
#   output.table = TRUE
#   match.by.chr = FALSE
#   duplicate.matches = c("best")
#   merge.matches = FALSE
#   minimum.match.length = 100
#   minimum.match.identity = 0.75
#   minimum.match.coverage = 0.75
#   add.flanks = 0
#   genome.search.string = NULL
#   threads = 8
#   memory = 32
#   overwrite = FALSE
#   quiet = FALSE
#   blast.path = conda.path


  #Same adds to blast path
  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

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
    genome.files = file.names[grep(".fai$", file.names, invert = T)]
    if (is.null(genome.search.string) != T){
      genome.files = file.names[grep(genome.search.string, file.names)]
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
    sample = gsub(pattern = ".fa$|.fna|.gz", replacement = "", x = genome.files[i])
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
      data.table::setnames(match.data, headers)

      #Matches need to be greater than 12
      filt.data = match.data[match.data$matches > minimum.match.length,]
      #Remove identical rows
      filt.data = filt.data[duplicated(filt.data) != T,]
      #Percent identitiy must match 50% or greater
      filt.data = filt.data[filt.data$pident >= minimum.match.identity,]

      #Make sure the hit is greater than 50% of the reference length
      filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$qLen),]

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
      }#end if

      #combines tables
      temp.filt = filt.data[!filt.data$tName %in% dup.names,]
      final.table = rbind(temp.filt, save.dup)
      final.table = final.table[order(final.table$tName)]

      #Extracts the genomic data using the final.table coordinates
      Rsamtools::indexFa(species.genome.path)
      fa = Rsamtools::FaFile(species.genome.path)
      gr = as(GenomicRanges::seqinfo(fa), "GRanges")

      #adds genome columns
      final.table[, gStart:=as.numeric(0)]
      final.table[, gEnd:=as.numeric(0)]
      header.data = colnames(final.table)

      for (j in 1:nrow(final.table)){
        #subsets data
        sub.match = final.table[j,]

        #Gets starts and stops
        gen.start = min(c(sub.match$qStart, sub.match$qEnd)) - add.flanks
        gen.end  = max(c(sub.match$qEnd, sub.match$qStart)) + add.flanks
        if (gen.start <= 0){ gen.start = as.numeric(1) }
        if (gen.end > sub.match$qLen) { gen.end = sub.match$qLen }
        if (gen.end-gen.start >= max(sub.match$tLen) * 100){ next }

        ## Save in table here
        data.table::set(final.table, i =  as.integer(j),
                        j = match("gStart", header.data), value = gen.start )

        data.table::set(final.table, i =  as.integer(j),
                        j = match("gEnd", header.data), value = gen.end )

      }#end j loop

      #removes skipped nas above
      final.table = final.table[is.na(final.table$gStart) != T,]

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
                    col.names = bed.headers,
                    quote = F,
                    sep = "\t")
      }#end if

      system(paste0("rm ", output.name, "/target_nucl-blast_db*"))

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
      gr = as(GenomicRanges::seqinfo(fa), "GRanges")

      #gets ranges of stuff and obtains sequences
      final.table = read.table(input.file, header = bed.headers)


      if (name.bed.names == TRUE){
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
        names(target.seqs) = paste0(gr.ranges$chrom, "_", gr.ranges$start, "-",gr.ranges$end, "_|_", sample)
      }#end else

    }#end bed

    #Writes the final loci
    final.loci = as.list(as.character(target.seqs))
    PhyloCap::writeFasta(sequences = final.loci, names = names(final.loci),
               paste0(species.dir, "/", sample, "_target-matches.fa"),
               nbchar = 1000000, as.string = T)

    print(paste0(sample, " Matching to the genome complete. ", length(final.loci), " targets extracted!"))

    if (length(zipped.up) == 1){ system(paste0("rm ", species.dir, "/", sample, "_genome.fa") ) }

  }# end i loop


} #End function


#END SCRIPT
