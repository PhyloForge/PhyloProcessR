#' @title removeContamination
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param input.reads path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.dir the new directory to save the adaptor trimmed sequences
#'
#' @param decontamination.path directory of genomes contaminants to scan samples
#'
#' @param mode "Sample" to run on a single sample or "Directory" to run on a directory of samples
#'
#' @param bbmap.path system path to bbmap in case it can't be found
#'
#' @param samtools.path system path to samtools in case it can't be found
#'
#' @param bwa.path system path to bwa in case it can't be found
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

removeContamination = function(input.reads = "adaptor-removed-reads",
                               output.dir = "decontaminated-reads",
                               decontamination.path = NULL,
                               mode = c("sample", "directory"),
                               map.match = 1,
                               read.mapper = "bwa",
                               bbmap.path = "bbsplit.sh",
                               samtools.path = "samtools",
                               bwa.path = "bwa",
                               threads = 1,
                               mem = 1,
                               resume = TRUE,
                               overwrite = FALSE,
                               quiet = TRUE) {

  #Debegging
  # decontamination.path = "/Users/chutter/Dropbox/Research/0_Github/Contamination_Genomes"
  # bbmap.path = "/usr/local/bin"
  # samtools.path = "/usr/local/bin/samtools"
  # bwa.path = "/usr/local/bin/bwa"
  # input.reads = "adaptor-removed-reads"
  # file.rename = "file_rename.csv"
  # output.dir = "decontaminated-reads"
  # mode = "directory"
  # threads = 4
  # mem = 8
  # resume = FALSE
  # overwrite = TRUE
  # quiet = TRUE
  # map.match = 0.99
  # read.mapper = "bwa"

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (is.null(decontamination.path) == TRUE){ stop("Please provide decontamination genomes / sequences.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.dir) == F){ dir.create(output.dir) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.dir))
      dir.create(output.dir)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Read in sample data
  reads = list.files(input.reads, recursive = T, full.names = T)
  sample.names = gsub(".*/", "", reads)
  sample.names = gsub("_R1_.*|_R2_.*|_READ1_.*|_READ2_.*", "", sample.names)
  sample.data = data.frame(File = sample.names, Sample = sample.names)

  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             rawReads =as.character(),
                             Task = as.character(),
                             Program = as.character(),
                             startPairs = as.numeric(),
                             removePairs = as.numeric())

  for (i in 1:nrow(sample.data)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    #Finds all files for this given sample and turns into a cat string
    sample.reads = reads[grep(pattern = paste0(sample.data$File[i], "_"), x = reads)]
    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = paste0(sample.data$Sample[i], "_"), x = reads)] }
    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.data$Sample[i], " does not have any reads present for files ",
           sample.data$File[i], " from the input spreadsheet. ")
    } #end if statement

    #################################################
    ### Part B: Create directories and move files
    #################################################
    #Create sample directory
    out.path = paste0(output.dir, "/", sample.data$Sample[i])
    dir.create(out.path)
    #Creates sample report directory
    report.path = paste0("logs/", sample.data$Sample[i])
    if (dir.exists(report.path) == FALSE) { dir.create(report.path) }

    #################################################
    ### Part C: Runs fastp
    #################################################
    #sets up output reads
    inread.1 = paste0(input.reads, "/", sample.data$Sample[i], "/", sample.data$Sample[i], "_READ1_L001.fastq.gz")
    inread.2 = paste0(input.reads, "/", sample.data$Sample[i], "/", sample.data$Sample[i], "_READ2_L001.fastq.gz")

    outread.1 = paste0(out.path, "/", sample.data$Sample[i], "_READ1_L001.fastq.gz")
    outread.2 = paste0(out.path, "/", sample.data$Sample[i], "_READ2_L001.fastq.gz")

    if (read.mapper == "bbsplit"){
      #Next runs bbsplit to remove other sources of contamination from other organisms
      system(paste0(bbmap.path, "/bbsplit.sh -Xmx", mem,"g in1=", inread.1, " in2=", inread.2,
                    " ref=", decontamination.path, " minid=", map.match,
                    " outu1=", outread.1, " outu2=", outread.2), ignore.stderr = T)
      system(paste0("rm -r ref"))
    }#end if

    if (read.mapper == "bwa"){
      #Create combined and indexed reference
      if (dir.exists("ref-index") == FALSE){
        #Create combined reference
        reference.list = list.files(decontamination.path)
        combined.ref = Biostrings::DNAStringSet()
        for (j in 1:length(reference.list)){
          ref.seq = Rsamtools::scanFa(Rsamtools::FaFile(paste0(decontamination.path, "/", reference.list[j])))   # loads up fasta file
          ref.seq = ref.seq[Biostrings::width(ref.seq) >= 100]
          names(ref.seq) = paste0(gsub(".fa$", "", reference.list[j]), "_seq-", rep(1:length(ref.seq)))
          combined.ref = append(combined.ref, ref.seq)
        }
        #Save reference
        dir.create("ref-index")
        ref.save = as.list(as.character(combined.ref))
        writeFasta(sequences = ref.save,
                   names = names(ref.save),
                   file.out = "ref-index/reference.fa")

        #To do: check if already exists
        #Indexes the reference
        system(paste0("bwa index -p ref-index/reference ref-index/reference.fa"))
      }#end dir ecists

      system(paste0(bwa.path, " mem -M -E -0 -k 100 -w 4 -L 100",
                    " -t ", threads, " ref-index/reference ",
                    inread.1, " ", inread.2, " | ", samtools.path ," sort -@ ", threads,
                    " -o ", out.path, "/decontam-all.bam"))

      #Need to figure out how to export reads back to normal
      system(paste0(samtools.path, " view -b -f 4 ", out.path, "/decontam-all.bam > ",
                    out.path, "/decontam-unmapped.bam"))
      system(paste0(samtools.path, " view -b -F 4 ", out.path, "/decontam-all.bam > ",
                    out.path, "/decontam-mapped.bam"))

      system(paste0(samtools.path, " fastq -@ ", threads, " ",
                    out.path, "/decontam-unmapped.bam -1 ", outread.1, " -2 ", outread.2))

      #Gets match statistics
      system(paste0(samtools.path, " index ", out.path, "/decontam-mapped.bam"))

      table.dat = (system(paste0(samtools.path, " idxstats ",
                                 out.path, "/decontam-mapped.bam | cut -f 1,3"), intern = T))
      table.dat = data.frame(Organism = gsub("\t.*", "", table.dat),
                             Count = as.numeric(gsub(".*\t", "", table.dat)) )
      table.dat$Organism = gsub("_.*", "", table.dat$Organism)

      contam.data = aggregate(table.dat$Count, FUN = sum, by = list(table.dat$Organism))
      contam.data = contam.data[-1,]
      colnames(contam.data) = c("Contaminant", "ReadPairs")

      system(paste0("rm ", out.path, "/decontam-*"))

      write.csv(contam.data, file = paste0("logs/", sample.data$Sample[i], "/contamination_summary.csv"), row.names = FALSE)

    }#end if

    #Gathers stats on initial data
    start.reads = as.numeric(system(paste0("zcat < ", inread.1, " | echo $((`wc -l`/4))"), intern = T))
    end.reads = as.numeric(system(paste0("zcat < ", outread.1, " | echo $((`wc -l`/4))"), intern = T))

    temp.remove = data.frame(Sample = sample.data$Sample[i],
                             rawReads = sample.data$File[i],
                             Task = "decontamination",
                             Program = read.mapper,
                             startPairs = start.reads,
                             removePairs = start.reads-end.reads,
                             endPairs = end.reads)

    summary.data = rbind(summary.data, temp.remove)

    print(paste0(sample.data$Sample[i], " Completed decontamination removal!"))
  }#end sample i loop

  write.csv(summary.data, file = paste0("logs/decontamination_summary.csv"), row.names = FALSE)

}#end function

