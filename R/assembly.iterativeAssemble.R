#' @title iterativeAssemble
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
iterativeAssemble = function(input.reads = NULL,
                             reference = NULL,
                             mapper = c("bwa", "hisat2"),
                             iterations = 5,
                             memory = 1,
                             threads = 1,
                             spades.path = NULL,
                             hisat2.path = NULL,
                             bwa.path = NULL,
                             cap3.path = NULL,
                             overwrite = FALSE,
                             quiet = TRUE) {

  #Debug
  setwd("/Volumes/LaCie/Mantellidae")
  input.reads = "/Volumes/LaCie/Mantellidae/Wakea_madinika_2001F54"
  reference = "/Volumes/LaCie/Ultimate_FrogCap/Final_Files/FINAL_marker-seqs_Mar14-2023.fa"
  iterations = 5
  spades.path = "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  samtools.path = "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  hisat2.path = "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  blast.path = "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  bwa.path = "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  cap3.path = "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  mapper = "bbmap"
  quiet = FALSE
  threads = 6

  if (is.null(blast.path) == FALSE){
    b.string = unlist(strsplit(blast.path, ""))
    if (b.string[length(b.string)] != "/") {
      blast.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { blast.path = "" }

  if (is.null(samtools.path) == FALSE){
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { samtools.path = "" }

  if (is.null(spades.path) == FALSE){
    b.string = unlist(strsplit(spades.path, ""))
    if (b.string[length(b.string)] != "/") {
      spades.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { spades.path = "" }

  #Same adds to bbmap path
  if (is.null(hisat2.path) == FALSE){
    b.string = unlist(strsplit(hisat2.path, ""))
    if (b.string[length(b.string)] != "/") {
      hisat2.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { hisat2.path = "" }

  #Same adds to bbmap path
  if (is.null(bwa.path) == FALSE){
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bwa.path = "" }

  #Same adds to bbmap path
  if (is.null(cap3.path) == FALSE){
    b.string = unlist(strsplit(cap3.path, ""))
    if (b.string[length(b.string)] != "/") {
      cap3.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { cap3.path = "" }

  #Quick checks
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (is.null(reference) == TRUE){ stop("Please provide a reference.") }

  #Writes reference to file if its not a file path
  if(file.exists("iterative_temp") == TRUE){ system(paste0("rm -r iterative_temp")) }
  if(file.exists("spades") == TRUE){ system(paste0("rm -r spades")) }
  dir.create("iterative_temp")

  #reference = Biostrings::readDNAStringSet(reference)

  #############################
  ## While loop start
  #############################
  for (i in 1:iterations){

    print(paste0("------------   iteration ", i, " begin ----------------------"))

    if (i == 1){ sensitivity = "-0.8" }
    if (i > 1){ sensitivity = "-0.2" }

    #Subsets the reads to the lane
    set.reads = paste0(input.reads, "/", list.files(input.reads))
    # set.read1 = set.reads[grep("READ1|R1", set.reads)]
    # if (length(input.reads) >= 2){ set.read2 = set.reads[grep("READ2|R2", set.reads)] }
    # if (length(input.reads) >= 3){ set.read3 = set.reads[grep("READ3|MERGE|singletons", set.reads)] }

    #Runs bbmap if selected
    if (mapper == "hisat2"){

      #Creates reference
      dir.create("iterative_temp/index")
      system(paste0(hisat2.path, "hisat2-build -f ", reference, " ",
                    "iterative_temp/index/reference"))

      if (length(set.reads) >= 2){

        ## Here is running hisat2
        system(paste0(hisat2.path, "hisat2 -q -x iterative_temp/index/reference",
                      " -1 ", set.reads[1], " -2 ", set.reads[2],
                      " -S iterative_temp/mapped_reads.sam --mp 1,0 --sp 1,0 --score-min L,0.0,", sensitivity,
                      " --threads ", threads))

        system(paste0(samtools.path, "samtools view -b -F 4 iterative_temp/mapped_reads.sam > iterative_temp/mapped_only.sam"))
        system(paste0(samtools.path, "samtools sort -n iterative_temp/mapped_only.sam > iterative_temp/mapped_sort.sam"))

        system(paste0(samtools.path, "samtools fastq -@ ", threads, " iterative_temp/mapped_sort.sam",
                      " -1 iterative_temp/sample_READ1.fastq.gz ",
                      " -2 iterative_temp/sample_READ2.fastq.gz "))

        system(paste0("rm iterative_temp/mapped_only.sam iterative_temp/mapped_reads.sam iterative_temp/mapped_sort.sam"))

        system("rm -r iterative_temp/index")
        temp.read.path = paste0("iterative_temp/", list.files("iterative_temp"))

      }#end 2 read


      if (length(set.reads) == 3){

        ## Here is running hisat2
        system(paste0(hisat2.path, "hisat2 -q -x iterative_temp/index/reference",
                      " -1 ", set.reads[1], " -2 ", set.reads[2], " -U ", set.reads[3],
                      " -S iterative_temp/mapped_reads.sam --mp 1,0 --sp 1,0 --score-min L,0.0,", sensitivity,
                      " --threads ", threads))

        system(paste0(samtools.path, "samtools view -b -F 4 iterative_temp/mapped_reads.sam > iterative_temp/mapped_only.sam"))
        system(paste0(samtools.path, "samtools sort -n iterative_temp/mapped_only.sam > iterative_temp/mapped_sort.sam"))

        system(paste0(samtools.path, "samtools fastq -@ ", threads, " iterative_temp/mapped_sort.sam",
                      " -1 iterative_temp/sample_READ1.fastq.gz ",
                      " -2 iterative_temp/sample_READ2.fastq.gz ",
                      " -s iterative_temp/sample_READ3.fastq.gz "))

        system(paste0("rm iterative_temp/mapped_only.sam iterative_temp/mapped_reads.sam iterative_temp/mapped_sort.sam"))

        system("rm -r iterative_temp/index")
        if (file.exists("iterative_temp/reference.fa")) { system("rm -r iterative_temp/reference.fa") }
        temp.read.path = paste0("iterative_temp/", list.files("iterative_temp"))

      }#end set read 3

    }#end bbmap if

    #Runs spades
    spades.contigs = runSpades(read.paths = temp.read.path,
                               full.path.spades = spades.path,
                               mismatch.corrector = F,
                               isolate = T,
                               quiet = T,
                               read.contigs = T,
                               threads = threads,
                               memory = memory)

    spades.contigs = spades.contigs[Biostrings::width(spades.contigs) >= 100]

    #Saves the raw reads themselves
    if (length(spades.contigs) == 0){ return(paste0("Sample ", input.reads, " failed, no assemblied contigs.")) }

    #Writes contigs for cap3
    write.loci = as.list(as.character(spades.contigs))
    PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                         "iterative_temp/blast_contigs.fa", nbchar = 1000000, as.string = T)

    #headers
    headers = c("qName", "tName", "pident", "matches", "misMatches", "gapopen",
                "qStart", "qEnd", "tStart", "tEnd", "evalue", "bitscore", "qLen", "tLen", "gaps")

    system(paste0(blast.path, "makeblastdb -in ", reference, " -parse_seqids -dbtype nucl",
                  " -out iterative_temp/blast_db"), ignore.stdout = quiet, ignore.stderr = quiet)

    #Matches samples to loci
    system(paste0(blast.path, "blastn -task dc-megablast -db iterative_temp/blast_db",
                  " -query iterative_temp/blast_contigs.fa -out iterative_temp/blast_match.txt",
                  " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen gaps\" ",
                  " -num_threads ", threads))

    #Need to load in transcriptome for each species and take the matching transcripts to the database
    if (length(readLines("iterative_temp/blast_match.txt")) == 0) { return(paste0("Sample ", input.reads, " failed, no assemblied contigs.")) }
    match.data = read.table("iterative_temp/blast_match.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
    colnames(match.data) = headers

    contigs = Biostrings::readDNAStringSet("iterative_temp/blast_contigs.fa")
    system(paste0("rm iterative_temp/blast_match.txt iterative_temp/blast_db* iterative_temp/blast_contigs.fa"))
    system(paste0("rm iterative_temp/sample_READ*"))
    save.contigs = contigs[names(contigs) %in% match.data$qName]
    mean(Biostrings::width(save.contigs))
    #Writes contigs for cap3
    write.loci = as.list(as.character(save.contigs))
    PhyloCap::writeFasta(sequences = write.loci, names = names(write.loci),
                         "iterative_temp/reference.fa", nbchar = 1000000, as.string = T)

    reference = "iterative_temp/reference.fa"

    #Prints completion info
    print(paste0("iteration ", i, " complete!"))
    print(paste0("-------------------------------------------------------"))


  }#end iterations if



    system("rm -r iterative_temp")
    return(combined.contigs)
  }#end else
  ##########################
}#end function






