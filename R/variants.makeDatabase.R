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

variants.makeDatabase = function(bam.directory = NULL,
                                 output.database = "variant-database",
                                 reference.path = "ref-index",
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
  # #subreference.name = "rag1"
  # output.database = "variant-database"
  # threads = 4 #number of threads, 8-10 is a recommended amount
  # memory = 8
  # gatk4.path = "/Users/chutter/miniconda3/bin"
  # resume = FALSE
  # quiet = FALSE
  # overwrite = TRUE

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
  if (dir.exists(output.database) == F){ dir.create(output.database) }

  #Get loci names I guess
  probe.loci = Biostrings::readDNAStringSet(paste0(reference.path, "/reference.fa"))
  loci.names = names(probe.loci)

  #Get multifile databases together
  sample.names = list.dirs(bam.directory, recursive = F, full.names = F)
  all.files = list.files(bam.directory, recursive = T, full.names = T)
  vcf.files = all.files[grep("gatk4-haplotype-caller.g.vcf.gz$", all.files)]
  vcf.string = paste0("-V ", vcf.files, collapse = " ")

  # #Resumes file download
  # if (resume == TRUE){
  #   done.files = list.files(output.directory, full.names = T, recursive = T)
  #   done.files = done.files[grep("gatk4-haplotype-caller.g.vcf.gz$", done.files)]
  #   done.names = gsub("/gatk4-haplotype-caller.g.vcf.gz$", "", done.files)
  #   done.names = gsub(".*\\/", "", done.names)
  #   sample.names = sample.names[!sample.names %in% done.names]
  # }
  #
  # #Resumes file download
  # if (overwrite == FALSE){
  #   done.files = list.files(output.directory, full.names = T, recursive = T)
  #   done.files = done.files[grep("gatk4-haplotype-caller.g.vcf.gz$", done.files)]
  #   done.names = gsub("/gatk4-haplotype-caller.g.vcf.gz$", "", done.files)
  #   done.names = gsub(".*\\/", "", done.names)
  #   sample.names = sample.names[!sample.names %in% done.names]
  # }

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
  foreach(i=1:length(loci.names), .packages = c("foreach")) %dopar% {

    #Loops through each locus and does operations on them
    #for (i in 1:length(loci.names)){

    #Combine them into a single database
    if (dir.exists(loci.names[i]) == T) { system(paste0("rm -r ", loci.names[i])) }
    system(paste0(gatk4.path, "gatk --java-options '-Xmx", mem.cl, "G' GenomicsDBImport ", vcf.string,
                  " --genomicsdb-workspace-path ", loci.names[i], " --intervals ", loci.names[i]))

    #Combine them into a single database
    system(paste0(gatk4.path, "gatk --java-options '-Xmx", mem.cl, "G' GenotypeGVCFs",
                  " -R ", reference.path, "/reference.fa -V gendb://", loci.names[i],
                  " --use-new-qual-calculator true -O ", output.database, "/", loci.names[i], ".vcf"))

    system(paste0("rm -r ", loci.names[i]))

  }#end i loop

  stopCluster(cl)

}#end function


#########################
###### END SCRIPT
#########################


