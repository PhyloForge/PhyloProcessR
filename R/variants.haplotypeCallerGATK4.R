#' @title variants.haplotypeCallerGATK4
#'
#' @description Function for running the program spades to assemble short read sequencing data
#'
#' @param read.directory directory of processed reads
#'
#' @param output.directory save name for the output directory
#'
#' @param full.path.spades contigs are added into existing alignment if algorithm is "add"
#'
#' @param mismatch.corrector algorithm to use: "add" add sequences with "add.contigs"; "localpair" for local pair align. All others available
#'
#' @param kmer.values if a file name is provided, save.name will be used to save aligment to file as a fasta
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

variants.haplotypeCallerGATK4 = function(bam.directory = NULL,
                                         output.directory = "variant-discovery",
                                         reference.path = "ref-index",
                                         samtools.path = NULL,
                                         bwa.path = NULL,
                                         picard.path = NULL,
                                         gatk4.path = NULL,
                                         threads = 1,
                                         memory = 1,
                                         resume = TRUE,
                                         overwrite = TRUE,
                                         quiet = TRUE) {

  #Debugging
  #Home directoroies
  # library(doParallel)
  # work.dir = "/Volumes/Armored/Test/variant-calling" #Your main project directory
  # dir.create(work.dir)
  # setwd(work.dir)
  #
  # bam.directory = "/Volumes/Armored/Test/variant-calling/variant-discovery"
  # reference.path = "ref-index"
  # subreference.name = "rag1"
  # output.directory = "variant-discovery"
  # auto.readgroup = T #Keep it T unless it crashes.
  # threads = 4 #number of threads, 8-10 is a recommended amount
  # memory = 8
  # samtools.path = "/Users/chutter/miniconda3/bin"
  # bwa.path = "/usr/local/bin"
  # picard.path = "/Users/chutter/miniconda3/bin"
  # gatk4.path = "/Users/chutter/miniconda3/bin"
  # resume = FALSE
  # quiet = FALSE
  # overwrite = TRUE

  #Same adds to bbmap path
  require(doParallel)
  if (is.null(samtools.path) == FALSE){
    b.string = unlist(strsplit(samtools.path, ""))
    if (b.string[length(b.string)] != "/") {
      samtools.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { samtools.path = "" }

  #Same adds to bbmap path
  if (is.null(bwa.path) == FALSE){
    b.string = unlist(strsplit(bwa.path, ""))
    if (b.string[length(b.string)] != "/") {
      bwa.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { bwa.path = "" }

  #Same adds to bbmap path
  if (is.null(picard.path) == FALSE){
    b.string = unlist(strsplit(picard.path, ""))
    if (b.string[length(b.string)] != "/") {
      picard.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { picard.path = "" }

  #Same adds to bbmap path
  if (is.null(gatk4.path) == FALSE){
    b.string = unlist(strsplit(gatk4.path, ""))
    if (b.string[length(b.string)] != "/") {
      gatk4.path = paste0(append(b.string, "/"), collapse = "")
    }#end if
  } else { gatk4.path = "" }

  #Quick checks
  if (is.null(bam.directory) == TRUE){ stop("Please provide the bam directory.") }
  if (file.exists(bam.directory) == F){ stop("BAM folder not found.") }

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Read in sample data
  bam.files = list.files(bam.directory, recursive = T, full.names = T)
  bam.files = bam.files[grep("final-mapped-all.bam$", bam.files)]
  sample.names = list.dirs(bam.directory, recursive = F, full.names = F)

  #Resumes file download
   if (resume == TRUE){
    done.files = list.files(output.directory, full.names = T, recursive = T)
    done.files = done.files[grep("gatk4-haplotype-caller.g.vcf.gz$", done.files)]
    done.names = gsub("/gatk4-haplotype-caller.g.vcf.gz$", "", done.files)
    done.names = gsub(".*\\/", "", done.names)
    sample.names = sample.names[!sample.names %in% done.names]
   }

  #Resumes file download
  if (overwrite == FALSE){
    done.files = list.files(output.directory, full.names = T, recursive = T)
    done.files = done.files[grep("gatk4-haplotype-caller.g.vcf.gz$", done.files)]
    done.names = gsub("/gatk4-haplotype-caller.g.vcf.gz$", "", done.files)
    done.names = gsub(".*\\/", "", done.names)
    sample.names = sample.names[!sample.names %in% done.names]
  }


  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################
  #Sets up multiprocessing
  cl = makeCluster(threads)
  registerDoParallel(cl)
  mem.cl = floor(memory/threads)

  #Loops through each locus and does operations on them
  foreach(i=1:length(sample.names), .packages = c("foreach", "ShortRead")) %dopar% {
    #Runs through each sample
    #for (i in 1:length(sample.names)) {
    #################################################
    ### Part A: prepare for loading and checks
    #################################################
    sample.dir = paste0(output.directory, "/", sample.names[i])
    if (file.exists(sample.dir) == FALSE) { dir.create(sample.dir) }

    #Gets the reads for the sample
    sample.bams = bam.files[grep(pattern = paste0(sample.names[i], "/"), x = bam.files)]
    #Checks the Sample column in case already renamed
    if (length(sample.bams) == 0){ sample.bams = sample.bams[grep(pattern = sample.names[i], x = sample.bams)] }

    #Returns an error if reads are not found
    if (length(sample.bams) == 0 ){
      stop(sample.names[i], " does not have any reads present for files ")
    } #end if statement

    #CReates new directory
    report.path = paste0("logs/", sample.names[i])
    if (file.exists(report.path) == FALSE) { dir.create(report.path) }

    #Sets up merging of bams from different lanes
    input.string = paste0("I=", sample.bams, collapse = " ")
    if (length(sample.bams) != 1){

      dir.create(paste0(sample.dir, "/Lane_Merge"))

      #Next combine .bam files together!
      system(paste0(picard.path, "picard -Xmx", mem.cl, "G MergeSamFiles",
                    " ", input.string, " O=", sample.dir, "/Lane_Merge/final-mapped-merge.bam",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

      #Sort by coordinate for input into MarkDuplicates
      system(paste0(picard.path, "picard -Xmx", mem.cl, "G SortSam",
                    " INPUT=", sample.dir, "/Lane_Merge/final-mapped-merge.bam",
                    " OUTPUT=", sample.dir, "/Lane_Merge/final-mapped-sort.bam",
                    " CREATE_INDEX=true SORT_ORDER=coordinate",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

      #Marks duplicate reads
      system(paste0(picard.path, "picard -Xmx", mem.cl, "G MarkDuplicates",
                    " INPUT=", sample.dir, "/Lane_Merge/final-mapped-sort.bam",
                    " OUTPUT=", sample.dir, "/Lane_Merge/final-mapped-dup.bam",
                    " CREATE_INDEX=true METRICS_FILE=logs/", sample.names[i], "/duplicate_metrics.txt",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

      #Sorts and stuff
      system(paste0(picard.path, "picard -Xmx", mem.cl, "G SortSam",
                    " INPUT=", sample.dir, "/Lane_Merge/final-mapped-dup.bam",
                    " OUTPUT=/dev/stdout SORT_ORDER=coordinate",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true | ",
                    picard.path, "picard -Xmx", mem.cl, "G SetNmAndUqTags",
                    " INPUT=/dev/stdin OUTPUT=", sample.dir, "/Lane_Merge/final-mapped-all.bam",
                    " CREATE_INDEX=true R=",reference.path, "/reference.fa",
                    " USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"))

      input.bam = paste0(sample.dir, "/Lane_Merge/final-mapped-all.bam")

      #Delete old files to make more space
      system(paste0("rm ", sample.dir, "/Lane_Merge/final-mapped-dup.bam"))
      system(paste0("rm ", sample.dir, "/Lane_Merge/final-mapped-sort.bam"))
      system(paste0("rm ", sample.dir, "/Lane_Merge/final-mapped-merge.bam"))

    } else {
      #lane 1 if thats all there is
      input.bam = paste0(sample.dir, "/Lane_1/final-mapped-all.bam")

    }# end else

    #Starts to finally look for Haplotypes! *here
    system(paste0(gatk4.path, "gatk --java-options '-Xmx", mem.cl, "G' HaplotypeCaller",
                  " -R ", reference.path, "/reference.fa -O ", sample.dir, "/gatk4-haplotype-caller.g.vcf.gz",
                  " -I ", input.bam,
                  " -ERC GVCF --max-alternate-alleles 3",
                  " -bamout ", sample.dir, "/gatk4-haplotype-caller.bam"))

    print(paste0(sample.names[i], " completed GATK4 haplotype caller!"))

  }# end i loop

  stopCluster(cl)


}#end function

# END SCRIPT
