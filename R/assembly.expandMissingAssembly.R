#' @title expandMissingAssembly
#'
#' @description Function for removing adaptor sequences from raw Illumina sequence data using the program fastp
#'
#' @param input.reads path to a folder of raw reads in fastq format.
#'
#' @param reference a csv file with a "File" and "Sample" columns, where "File" is the file name and "Sample" is the desired renamed file
#'
#' @param output.name the new directory to save the adaptor trimmed sequences
#'
#' @param mapper "Sample" to run on a single sample or "Directory" to run on a directory of samples
#'
#' @param min.iterations system path to fastp in case it can't be found
#'
#' @param max.iterations system path to fastp in case it can't be found
#'
#' @param min.length system path to fastp in case it can't be found
#'
#' @param max.length system path to fastp in case it can't be found
#'
#' @param min.ref.id system path to fastp in case it can't be found
#'
#' @param spades.path system path to fastp in case it can't be found
#'
#' @param bbmap.path system path to fastp in case it can't be found
#'
#' @param cap3.path system path to fastp in case it can't be found
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

#Iteratively assembles to reference
expandMissingAssembly = function(assembly.directory = NULL,
                                 read.directory = NULL,
                                 reference = NULL,
                                 output.directory = "expand-missing-assembly",
                                 mapper = c("bwa", "hisat2"),
                                 memory = 1,
                                 threads = 1,
                                 spades.path = NULL,
                                 hisat2.path = NULL,
                                 bwa.path = NULL,
                                 cap3.path = NULL,
                                 overwrite = FALSE,
                                 quiet = TRUE) {

  #Debug
  library(PhyloProcessR)
  setwd("/Volumes/LaCie/Mantellidae")
  assembly.directory = "/Volumes/LaCie/Mantellidae/draft-assemblies"
  output.directory = "expand-missing-assembly"
  reference = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
  read.directory = "/Volumes/LaCie/Mantellidae/reads"

  spades.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  samtools.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  hisat2.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  blast.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  bwa.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  fastp.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  cdhit.path = "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  mapper = "bbmap"
  quiet = TRUE
  overwrite = FALSE
  threads = 8
  memory = 1024
  min.match.percent = 60
  min.match.length = 100
  min.match.coverage = 35

  if (is.null(fastp.path) == FALSE) {
    b.string <- unlist(strsplit(fastp.path, ""))
    if (b.string[length(b.string)] != "/") {
      fastp.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    fastp.path <- ""
  }

  if (is.null(blast.path) == FALSE) {
    b.string <- unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    blast.path <- ""
  }

  if (is.null(samtools.path) == FALSE) {
    b.string <- unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    samtools.path <- ""
  }

  if (is.null(spades.path) == FALSE) {
    b.string <- unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    spades.path <- ""
  }

  # Same adds to bbmap path
  if (is.null(hisat2.path) == FALSE) {
    b.string <- unlist(strsplit(hisat2.path, ""))
    if (b.string[length(b.string)] != "/") {
      hisat2.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    hisat2.path <- ""
  }

  # Same adds to bbmap path
  if (is.null(bwa.path) == FALSE) {
    b.string <- unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    bwa.path <- ""
  }

  # Same adds to bbmap path
  if (is.null(cdhit.path) == FALSE) {
    b.string <- unlist(strsplit(cdhit.path, ""))
    if (b.string[length(b.string)] != "/") {
      cdhit.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    cdhit.path <- ""
  }

  # Quick checks
  if (is.null(assembly.directory) == TRUE) {
    stop("Please provide input reads.")
  }
  if (is.null(reference) == TRUE) {
    stop("Please provide a reference.")
  }

  # #Writes reference to file if its not a file path
  # if(file.exists("iterative_temp") == TRUE){ system(paste0("rm -r iterative_temp")) }
  # if(file.exists("spades") == TRUE){ system(paste0("rm -r spades")) }
  # dir.create("iterative_temp")

  reference.seqs = Biostrings::readDNAStringSet(reference)

  file.names = list.files(assembly.directory)
  read.sets = list.files(read.directory)

  dir.create(output.directory)

  #############################
  ## Target matching loop start
  #############################

  #headers
  headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
              "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

  save.samples = c()
  for (i in seq_along(file.names)){

    #Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])
    species.dir = paste0(output.directory, "/", sample)

    # Creates species directory if none exists
    if (file.exists(species.dir) == FALSE) {
      dir.create(species.dir)
    }

    # Make blast database for the probe loci
    system(paste0(
      blast.path, "makeblastdb -in ", reference,
      " -parse_seqids -dbtype nucl -out ", species.dir, "/nucl-blast_db"
    ), ignore.stdout = quiet)

    # Matches samples to loci
    system(paste0(
      blast.path, "blastn -task dc-megablast -db ", species.dir, "/nucl-blast_db -evalue 0.001",
      " -query ", assembly.directory, "/", file.names[i], " -out ", species.dir, "/target-blast-match.txt",
      " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
      " -num_threads ", threads
    ))

    # Need to load in transcriptome for each species and take the matching transcripts to the database
    system(paste0("rm ", species.dir, "/nucl-blast_db*"))

    # Loads in match data
    match.data = data.table::fread(paste0(species.dir, "/target-blast-match.txt"),
      sep = "\t", header = FALSE, stringsAsFactors = FALSE
    )
    data.table::setnames(match.data, headers)

    #Matches need to be greater than 12
    filt.data = match.data[match.data$matches > min.match.length,]
    #Percent identitiy must match 50% or greater
    filt.data = filt.data[filt.data$pident >= min.match.percent,]

    #Make sure the hit is greater than 50% of the reference length
    filt.data = filt.data[filt.data$matches >= ((min.match.coverage / 100) * filt.data$tLen), ]

    if (nrow(filt.data) == 0) {
      print(paste0(sample, " had no matches. Skipping"))
      next
    }

    #Sorting: exon name, contig name, bitscore higher first, evalue
    data.table::setorder(filt.data, qName, tName, -pident, -bitscore, evalue)

    sample.contigs = Biostrings::readDNAStringSet(paste0(assembly.directory, "/", file.names[i]))
    match.contigs = sample.contigs[names(sample.contigs) %in% filt.data$qName]

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(match.contigs))
    writeFasta(sequences = final.loci, names = names(final.loci),
              paste0(species.dir, "/", sample, "_matching-contigs.fa"), nbchar = 1000000, as.string = TRUE)


    system(paste0("rm ", species.dir, "/target-blast-match.txt"))

    write.table(filt.data, file = paste0(species.dir, "/filtered-blast-match.txt"),
                row.names = FALSE, quote = FALSE, sep = "\t")


    # miss.reference = reference.seqs[!names(reference.seqs) %in% filt.data$qName]
    #
    # write.loci = as.list(as.character(miss.reference))
    # PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
    #                      paste0(species.dir, "/missing_ref.fa"), nbchar = 1000000, as.string = T)
    #
    # reference = paste0(species.dir, "/missing_ref.fa")

    final.match = c()
    for (j in seq_len(nrow(filt.data))){
      temp.match = match.contigs[names(match.contigs) %in% filt.data$qName[j]]
      names(temp.match) = filt.data$tName[j]
      final.match = append(final.match, temp.match)
    }#end j loop

    save.samples = append(save.samples, final.match)

  }#End target matching loop

  dup.names = unique(names(save.samples)[duplicated(names(save.samples)) == TRUE])
  non.dup = save.samples[!names(save.samples) %in% dup.names]

  final.save = c()
  for (j in seq_along(dup.names)){
    dup.contigs = save.samples[names(save.samples) %in% dup.names[j]]
    save.long = dup.contigs[Biostrings::width(dup.contigs) == max(Biostrings::width(dup.contigs))][1]
    final.save = append(final.save, save.long)
  }

  final.save = append(final.save, non.dup)

  #Finds probes that match to two or more contigs
  final.loci = as.list(as.character(final.save))
  writeFasta(sequences = final.loci, names = names(final.loci),
            paste0(output.directory, "/unique_matches.fa"), nbchar = 1000000, as.string = TRUE)

  #####################################################################################
  #### Read mapping and assembly step
  #####################################################################################

  ### Approach 1
  # 1. take assembled matching contigs for each sample
    # a. don't modify anything, rename contigs whether duplicated or not
  # 2. map reads to assembled contigs
  # 3. extract reads that don't map, assemble those

  ### Approach 2.
  # 1. take assembled matching contigs for each sample
  # a. don't modify anything, rename contigs whether duplicated or not
  # 2. combine all contigs for all samples
  # 3. remove duplicate names, choosing longest
  # 4. map reads to unique sample database
  # 5. assemble

  dir.create("expanded-assemblies")

  for (i in seq_along(file.names)){

    #Sets up working directories for each species
    sample = gsub(pattern = ".fa$", replacement = "", x = file.names[i])
    species.dir <- paste0(output.directory, "/", sample)

    # Checks if this has been done already
    if (overwrite == FALSE) {
      if (file.exists(paste0("expanded-assemblies/", sample, ".fa")) == TRUE) {
        print(paste0(sample, " already finished, skipping. Set overwrite to T to overwrite."))
        next
      }
    } # end

    # Gathers reads
    input.reads = read.sets[read.sets == sample]
    set.reads = list.files(paste0(read.directory, "/", input.reads), full.names = TRUE)

    # found data
    found.data = read.table(paste0(species.dir, "/filtered-blast-match.txt"),
      sep = "\t", header = TRUE, stringsAsFactors = FALSE
    )

    # found.data1 = read.table(paste0(species.dir, "/found-missing-blast-match.txt"), sep = "\t", header = T, stringsAsFactors = FALSE)
    # found.data = rbind(found.data1, found.data2)

    sample.save = final.save[!names(final.save) %in% found.data$tName]

    # Finds probes that match to two or more contigs
    final.loci = as.list(as.character(sample.save))
    writeFasta(
      sequences = final.loci, names = names(final.loci),
      paste0(species.dir, "/missing_ref.fa"), nbchar = 1000000, as.string = T
    )

    missing.ref = paste0(species.dir, "/missing_ref.fa")

    # missing.ref = "/Volumes/LaCie/Mantellidae/find-missing-asssemble/Blommersia_grandisonae_CRH-792/Blommersia_grandisonae_CRH-792_matching-contigs.fa"

    #Creates reference
    dir.create(paste0(species.dir, "/index"))
    system(paste0(hisat2.path, "hisat2-build -f ", missing.ref, " ",
                  species.dir, "/index/reference"))

    # system(paste0(hisat2.path, "hisat2 -q -x ", species.dir, "/index/reference",
    #               " -1 ", set.reads[1], " -2 ", set.reads[2],
    #               " -S ", species.dir, "/mapped_reads.sam ",
    #               " --threads ", threads))

    # sample.reads <- unique(gsub("_1.f.*|_2.f.*|_3.f.*|-1.f.*|-2.f.*|-3.f.*|_R1_.*|_R2_.*|_R3_.*|_READ1_.*|_READ2_.*|_READ3_.*|_R1.f.*|_R2.f.*|_R3.f.*|-R1.f.*|-R2.f.*|-R3.f.*|_READ1.f.*|_READ2.f.*|_READ3.f.*|-READ1.f.*|-READ2.f.*|-READ3.f.*|_singleton.*|-singleton.*|READ-singleton.*|READ_singleton.*|_READ-singleton.*|-READ_singleton.*|-READ-singleton.*|_READ_singleton.*", "", sample.reads))

    # Here is running hisat2
    system(paste0(hisat2.path, "hisat2 -q -x ", species.dir, "/index/reference",
                  " -1 ", set.reads[1], " -2 ", set.reads[2],
                  " -S ", species.dir, "/mapped_reads.sam --mp 1,0 --sp 1,0 --score-min L,0.0,-0.3",
                  " --threads ", threads))

    system(paste0(
      samtools.path, "samtools view -@ ", threads, " -b -F 4 ",
      species.dir, "/mapped_reads.sam > ", species.dir, "/mapped_all.sam"
    ))
    system(paste0(
      samtools.path, "samtools view -@ ", threads, " -b -f 4 -F 264 ",
      species.dir, "/mapped_reads.sam > ", species.dir, "/mapped1.sam"
    ))
    system(paste0(
      samtools.path, "samtools view -@ ", threads, " -b -f 8 -F 260 ",
      species.dir, "/mapped_reads.sam > ", species.dir, "/mapped2.sam"
    ))
    system(paste0(
      samtools.path, "samtools merge ", species.dir, "/mapped_combined.sam ",
      species.dir, "/mapped_all.sam ", species.dir, "/mapped1.sam ", species.dir, "/mapped2.sam"
    ))
    system(paste0(
      samtools.path, "samtools sort -n ",
      species.dir, "/mapped_combined.sam > ", species.dir, "/mapped_sort.sam"
    ))

    # system(paste0(samtools.path, "samtools flagstat ", species.dir, "/mapped_sort.sam"))
    #
    # system(paste0(samtools.path, "samtools view -@ ", threads, " -b -f 12 -F 256 ", species.dir, "/mapped_reads.sam > ", species.dir, "/unmapped_only.sam"))
    # system(paste0(samtools.path, "samtools sort -n ", species.dir, "/unmapped_only.sam > ", species.dir, "/unmapped_sort.sam"))
    #
    dir.create(paste0(species.dir, "/temp_reads"))
    dir.create(paste0(species.dir, "/temp_reads/sample"))
    system(paste0(samtools.path, "samtools fastq -@ ", threads, " ", species.dir, "/mapped_sort.sam",
                  " -1 ", species.dir, "/temp_reads/sample/sample_READ1.fastq.gz ",
                  " -2 ", species.dir, "/temp_reads/sample/sample_READ2.fastq.gz "),
           ignore.stdout = quiet, ignore.stderr = quiet)

    #Merge paired-end reads for de novo assembly
    mergePairedEndReads(input.reads = paste0(species.dir, "/temp_reads"),
                        output.directory =  paste0(species.dir, "/temp-merged-reads"),
                        fastp.path = fastp.path,
                        threads = threads,
                        mem = memory,
                        overwrite = overwrite,
                        quiet = FALSE)

    #Removes if nothing merged
    if (file.size(paste0(species.dir, "/temp-merged-reads/sample/sample_READ3.fastq.gz")) == 0){
      system(paste0("rm ", species.dir, "/temp-merged-reads/sample/sample_READ3.fastq.gz"))
    }

    system(paste0("rm ", species.dir, "/mapped_all.sam ",
                  species.dir, "/mapped_combined.sam ",
                  species.dir, "/mapped_reads.sam ",
                  species.dir, "/mapped_sort.sam ",
                  species.dir, "/mapped1.sam ",
                  species.dir, "/mapped2.sam"))

    temp.read.path = paste0(species.dir, "/temp-merged-reads/sample/", list.files(paste0(species.dir, "/temp-merged-reads/sample")))
    system(paste0("rm -rf ", species.dir, "/temp_reads"))

    #Runs spades
    spades.contigs = runSpades(read.paths = temp.read.path,
                              full.path.spades = spades.path,
                              mismatch.corrector = TRUE,
                              isolate = FALSE,
                              quiet = quiet,
                              read.contigs = TRUE,
                              threads = threads,
                              memory = memory)

    spades.contigs = spades.contigs[Biostrings::width(spades.contigs) >= 100]

    #Saves the raw reads themselves
    if (length(spades.contigs) == 0){ return(paste0("Sample ", input.reads, " failed, no assemblied contigs.")) }

    #Writes contigs for cap3
    write.loci = as.list(as.character(spades.contigs))
    PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                         paste0(species.dir, "/blast_contigs.fa"), nbchar = 1000000, as.string = T)

    #headers
    headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

    system(paste0(blast.path, "makeblastdb -in ", missing.ref, " -parse_seqids -dbtype nucl",
                  " -out ", species.dir, "/blast_db"), ignore.stdout = quiet, ignore.stderr = quiet)

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db ", species.dir, "/blast_db",
                  " -query ", species.dir, "/blast_contigs.fa -out ", species.dir, "/blast_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))

    #Need to load in transcriptome for each species and take the matching transcripts to the database
    if (length(readLines(paste0(species.dir, "/blast_match.txt"))) == 0) { return(paste0("Sample ", input.reads, " failed, no new assemblied contigs.")) }
    match.data = read.table(paste0(species.dir, "/blast_match.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    colnames(match.data) = headers

    sample.contigs = Biostrings::readDNAStringSet(paste0(species.dir, "/blast_contigs.fa"))
    system(paste0("rm ", species.dir, "/blast_match.txt ",  species.dir, "/blast_db* ",
                  species.dir, "/blast_contigs.fa"))
    system(paste0("rm -r ", species.dir, "/temp-merged-reads"))

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
    filt.data = filt.data[filt.data$matches >= ( (min.match.coverage/100) * filt.data$tLen),]

    match.contigs = sample.contigs[names(sample.contigs) %in% filt.data$qName]

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(match.contigs))
    writeFasta(sequences = final.loci, names = names(final.loci),
               paste0(species.dir, "/", sample, "_found-contigs.fa"), nbchar = 1000000, as.string = T)

    write.table(filt.data, file = paste0(species.dir, "/found-missing-blast-match.txt"),
                row.names = F, quote = F, sep = "\t")

    print(paste0(sample, " target matching complete. ", length(unique(found.data$tName)), " targets found in the first round. ",
                 length(match.contigs), " found in the rebuild missing round!"))


    #Saves final contigs
    found.contigs = Biostrings::readDNAStringSet(paste0(species.dir, "/", sample, "_matching-contigs.fa"))
    output.contigs = append(match.contigs, found.contigs)

    #Finds probes that match to two or more contigs
    final.loci = as.list(as.character(output.contigs))
    writeFasta(sequences = final.loci, names = names(final.loci),
               paste0("expanded-assemblies/", sample, ".fa"), nbchar = 1000000, as.string = T)

    #0.6 = 498
    #0.5 = 495
    #0.4 = 481
    #0.3 = 443
    #0.2 = 382

  }#end iterations if
  ##########################
}#end function






