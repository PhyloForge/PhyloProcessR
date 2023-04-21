#' @title baseRecalibration
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

baseRecalibration = function(haplotype.caller.directory = "haplotype-caller",
                            mapping.directory = "sample-mapping",
                            gatk4.path = NULL,
                            threads = 1,
                            memory = 1,
                            clean.up = TRUE,
                            overwrite = TRUE,
                            quiet = TRUE) {

  #Debugging
  #Home directoroies
  # library(PhyloCap)
  # setwd("/Volumes/LaCie/Mantellidae/data-analysis")
  # haplotype.caller.directory <- "variant-calling/haplotype-caller"
  # mapping.directory <- "variant-calling/sample-mapping"

  # gatk4.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"
  # samtools.path <- "/Users/chutter/Bioinformatics/anaconda3/envs/PhyloCap/bin"

  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE
  # clean.up = TRUE

  # Same adds to bbmap path
  require(foreach)

  # Same adds to bbmap path
  if (is.null(gatk4.path) == FALSE) {
    b.string <- unlist(strsplit(gatk4.path, ""))
    if (b.string[length(b.string)] != "/") {
      gatk4.path <- paste0(append(b.string, "/"), collapse = "")
    } # end if
  } else {
    gatk4.path <- ""
  }

  #Quick checks
  if (is.null(haplotype.caller.directory) == TRUE) {
    stop("Please provide the haplotype caller directory.")
  }
  if (file.exists(haplotype.caller.directory) == F) {
    stop("Haplotype caller directory not found.")
  }

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  #Get multifile databases together
  sample.names = list.dirs(haplotype.caller.directory, recursive = F, full.names = F)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(haplotype.caller.directory, full.names = T, recursive = T)
    done.files <- done.files[grep("gatk4-bqsr-haplotype-caller.g.vcf.gz$", done.files)]
    done.names <- gsub("/gatk4-bqsr-haplotype-caller.g.vcf.gz$", "", done.files)
    done.names <- gsub(".*\\/", "", done.names)
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0){ return("no samples available to analyze.") }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  # Sets up multiprocessing
  cl <- snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl <- floor(memory / threads)

  #Loops through each locus and does operations on them
  foreach::foreach(i=seq_along(sample.names), .packages = c("foreach")) %dopar% {
    # Loops through each locus and does operations on them
    # for (i in 1:length(loci.names)){

    reference.path <- paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " GenotypeGVCFs -R ", reference.path,
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-haplotype-caller.g.vcf.gz",
      " --use-new-qual-calculator true -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-genotype.vcf"
    ))

    # Selects only the SNPs from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-genotype.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf",
      " --select-type SNP"
    ))

    #Selects only the indels from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-genotype.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf",
      " --select-type INDEL"
    ))

    # Applies filters to SNPs
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf",
      " -filter \"QD<2.0\" --filter-name \"QD2\"",
      " -filter \"QUAL<100.0\" --filter-name \"QUAL30\"",
      " -filter \"SOR>3.0\" --filter-name \"SOR3\"",
      " -filter \"FS>60.0\" --filter-name \"FS60\"",
      " -filter \"MQ<50.0\" --filter-name \"MQ40\"",
      " -filter \"MQRankSum<-12.5\" --filter-name \"MQRankSum12.5\"",
      " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum8\""
    ))

    # Applies filters to indels
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf",
      " -filter \"QD<2.0\" --filter-name \"QD2\"",
      " -filter \"QUAL<100.0\" --filter-name \"QUAL30\"",
      " -filter \"FS>200.0\" --filter-name \"FS200\"",
      " -filter \"ReadPosRankSum<-20.0\" --filter-name \"ReadPosRankSum20\""
    ))

    # Combine them into a single VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " SortVcf",
      " -I ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf",
      " -I ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-combined.vcf"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " SelectVariants",
      " -V ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-combined.vcf",
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf",
      " --exclude-filtered TRUE"
    ))

    #Gathers correct read files
    lane.files = list.dirs(paste0(mapping.directory, "/", sample.names[i]))
    lane.files = lane.files[grep("Lane_", lane.files)]

    if (length(lane.files) != 1) {
      read.bam = paste0(mapping.directory, "/", sample.names[i], "/Lane_Merge")
    }
    if (length(lane.files) == 1) {
      read.bam = paste0(mapping.directory, "/", sample.names[i], "/Lane_1")
    }

    #Applies base recalibration to the original bam file
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " BaseRecalibrator",
      " -I ", read.bam, "/final-mapped-all.bam",
      " -R ", reference.path,
      " --known-sites ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf",
      " -O ", read.bam, "/recal_data.table"
    ))

    recal.bam <- paste0(read.bam, "/bqsr-mapped-all.bam")

    # Applies the recalibration and creates a new bam file
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " ApplyBQSR",
      " -I ", read.bam, "/final-mapped-all.bam",
      " -R ", reference.path,
      " --bqsr-recal-file ", read.bam, "/recal_data.table",
      " -O ", recal.bam
    ))

   # Starts to finally look for Haplotypes! *here
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", mem.cl, "G\"",
      " HaplotypeCaller",
      " -R ", reference.path,
      " -I ", recal.bam,
      " -O ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-haplotype-caller.g.vcf.gz",
      " -ERC GVCF",
      " -bamout ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-haplotype-caller.bam"
    ))

    if (clean.up == TRUE){

      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-combined.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-indels.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-filtered-snps.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-rem-filtered-combined.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-snps.vcf.idx"))
      system(paste0("rm ", haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-indels.vcf.idx"))

    }

    print(paste0(sample.names[i], " completed GATK4 base recalibration!"))

  }#end i loop

  snow::stopCluster(cl)

}#end function

# END SCRIPT
