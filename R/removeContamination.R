#' @title removeContamination
#'
#' @description Function for removing contamination from other organisms from adaptor trimmed Illumina sequence data using BWA
#'
#' @param input.reads path to a folder of adaptor trimmed reads in fastq format.
#'
#' @param output.directory the new directory to save the adaptor trimmed sequences
#'
#' @param decontamination.path directory of genomes contaminants to scan samples
#'
#' @param mode "Sample" to run on a single sample or "Directory" to run on a directory of samples
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
                               output.directory = "decontaminated-reads",
                               decontamination.path = NULL,
                               mode = c("sample", "directory"),
                               map.match = 1,
                               samtools.path = "samtools",
                               bwa.path = "bwa",
                               threads = 1,
                               mem = 1,
                               resume = TRUE,
                               overwrite = FALSE,
                               quiet = TRUE) {

  #Debug
  # setwd("/Users/chutter/Dropbox/Research/0_Github/Test-dataset")
  # read.directory = "read-processing/adaptor-removed-reads"
  # output.directory = "read-processing/decontaminated-reads"
  # decontamination.path = "/Users/chutter/Dropbox/Research/0_Github/Contamination_Genomes"
  # samtools.path = "/Users/chutter/miniconda3/bin/samtools"
  # bwa.path = "/usr/local/bin/bwa"
  # mode = "directory"
  # threads = 4
  # mem = 8
  # resume = TRUE
  # overwrite = FALSE
  # quiet = TRUE
  # map.match = 0.99

  #Quick checks
  options(stringsAsFactors = FALSE)
  if (is.null(input.reads) == TRUE){ stop("Please provide input reads.") }
  if (file.exists(input.reads) == F){ stop("Input reads not found.") }
  if (is.null(decontamination.path) == TRUE){ stop("Please provide decontamination genomes / sequences.") }

  #Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == F){ dir.create(output.directory) } else {
    if (overwrite == TRUE){
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  }#end else

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Read in sample data **** sample is run twice?!
  reads = list.files(input.reads, recursive = T, full.names = T)
  sample.names = list.files(input.reads, recursive = F, full.names = F)

  #Resumes file download
  if (resume == TRUE){
    done.files = list.files(output.directory)
    sample.names = sample.names[!sample.names %in% done.files]
  }

  if (length(sample.names) == 0){ stop("no samples remain to analyze.") }

  #Creates the summary log
  summary.data =  data.frame(Sample = as.character(),
                             Lane = as.character(),
                             Task = as.character(),
                             Program = as.character(),
                             startPairs = as.numeric(),
                             removePairs = as.numeric())

  #Runs through each sample
  for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.reads = reads[grep(pattern = paste0(sample.names[i], "_"), x = reads)]
    sample.reads = unique(gsub("_R1_.*|_R2_.*|_READ1_.*|_READ2_.*|_R1.fast.*|_R2.fast.*|_READ1.fast.*|_READ2.fast.*", "", sample.reads))

    #Checks the Sample column in case already renamed
    if (length(sample.reads) == 0){ sample.reads = reads[grep(pattern = paste0(sample.names[i], x = reads))] }
    #Returns an error if reads are not found
    if (length(sample.reads) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    #CReates new directory
    out.path = paste0(output.directory, "/", sample.names[i])
    report.path = paste0("logs/", sample.names[i])
    if (file.exists(out.path) == FALSE) { dir.create(out.path) }
    if (file.exists(report.path) == FALSE) { dir.create(report.path) }

    for (j in 1:length(sample.reads)){

      lane.reads = reads[grep(pattern = paste0(sample.reads[j], "_"), x = reads)]

      #Checks the Sample column in case already renamed
      if (length(lane.reads) == 0){ lane.reads = reads[grep(pattern = paste0(sample.reads[j], x = reads))] }
      #Returns an error if reads are not found
      if (length(lane.reads) == 0 ){
        stop(sample.reads[j], " does not have any reads present for files ")
      } #end if statement

      lane.save = gsub(".*/", "", lane.reads)
      lane.name = gsub(".*/", "", sample.reads[j])

      #################################################
      ### Part C: Runs fastp
      #################################################
      #sets up output reads
      outreads = paste0(out.path, "/", lane.save)

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
        system(paste0(bwa.path, " index -p ref-index/reference ref-index/reference.fa"))
      }#end dir ecists

      #### REMOVE BOTH READ PAIRS

      system(paste0(bwa.path, " mem -M -E -0 -k 100 -w 4 -L 100",
                    " -t ", threads, " ref-index/reference ",
                    lane.reads[1], " ", lane.reads[2], " | ", samtools.path ," sort -@ ", threads,
                    " -o ", out.path, "/decontam-all.bam"))

      #Need to figure out how to export reads back to normal
      system(paste0(samtools.path, " view -b -f 4 ", out.path, "/decontam-all.bam > ",
                    out.path, "/decontam-unmapped.bam"))
      system(paste0(samtools.path, " view -b -F 4 ", out.path, "/decontam-all.bam > ",
                    out.path, "/decontam-mapped.bam"))

      system(paste0(samtools.path, " fastq -@ ", threads, " ",
                    out.path, "/decontam-unmapped.bam -1 ", outreads[1], " -2 ", outreads[2]))

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

      #Gathers stats on initial data
      start.reads = as.numeric(system(paste0("zcat < ", lane.reads[1], " | echo $((`wc -l`/4))"), intern = T))
      end.reads = as.numeric(system(paste0("zcat < ", outreads[1], " | echo $((`wc -l`/4))"), intern = T))

      temp.remove = data.frame(Sample = sample.names[i],
                               Lane = gsub(".*_", "", lane.name),
                               Task = "decontamination",
                               Program = "bwa",
                               startPairs = start.reads,
                               removePairs = start.reads-end.reads,
                               endPairs = end.reads)

      summary.data = rbind(summary.data, temp.remove)
      write.csv(contam.data, file = paste0("logs/", sample.names[i], "/contamination-read-counts.csv"), row.names = FALSE)

      print(paste0(lane.name, " Completed decontamination removal!"))

    }#end sample j loop

    print(paste0(sample.names[i], " Completed decontamination removal!"))

  }#end sample i loop

  write.csv(summary.data, file = paste0("logs/removeContamination_summary.csv"), row.names = FALSE)

}

