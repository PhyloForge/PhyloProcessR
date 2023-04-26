#' @title genotypeSamples
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

variants.genotypeSamples = function(mapping.directory = "sample-mapping",
                                    haplotype.caller.directory = "haplotype-caller",
                                    output.directory = "sample-genotypes",
                                    filtering.thresholds = c("high", "medium", "low"),
                                    gatk4.path = NULL,
                                    threads = 1,
                                    memory = 1,
                                    overwrite = TRUE,
                                    quiet = TRUE) {

 #Debugging
  #Home directoroies
  # library(PhyloCap)
  # setwd("/Volumes/LaCie/Mantellidae")
  # variant.calling.directory = "variant-calling"
  # output.directory <- "sample-genotypes"
  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"

  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- TRUE

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

  # Quick checks
  if (is.null(mapping.directory) == TRUE) {
    stop("Please provide the sample mapping directory.")
  }
  if (file.exists(mapping.directory) == FALSE) {
    stop("Sample mapping directory not found.")
  }

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
    }
  } # end else

  # Get multifile databases together
  sample.names <- list.dirs(haplotype.caller.directory, recursive = F, full.names = F)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(output.directory, full.names = T, recursive = T)
    done.files <- done.files[grep("gatk4-final-genotypes.vcf$", done.files)]
    done.names <- gsub("/gatk4-final-genotypes.vcf$", "", done.files)
    done.names <- gsub(".*\\/", "", done.names)
    sample.names <- sample.names[!sample.names %in% done.names]
  }

  if (length(sample.names) == 0) {
    return("no samples available to analyze.")
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
  foreach::foreach(i = seq_along(sample.names), .packages = c("foreach")) %dopar% {
    # Loops through each locus and does operations on them
    # for (i in 1:length(loci.names)){

    dir.create(paste0(output.directory, "/", sample.names[i]))

    reference.path <- paste0(mapping.directory, "/", sample.names[i], "/index/reference.fa")

    if (file.exists(paste0(haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-haplotype-caller.g.vcf.gz")) == TRUE ) {
      haplocaller.file <- paste0(haplotype.caller.directory, "/", sample.names[i], "/gatk4-bqsr-haplotype-caller.g.vcf.gz")
    } else {
      haplocaller.file <- paste0(haplotype.caller.directory, "/", sample.names[i], "/gatk4-haplotype-caller.g.vcf.gz")
    }

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
      " GenotypeGVCFs -R ", reference.path,
      " -V ", haplocaller.file,
      " --use-new-qual-calculator true -O ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-genotypes.vcf"
    ))

    # Selects only the SNPs from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-genotypes.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-snps.vcf",
      " --select-type SNP"
    ))

    #Selects only the indels from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-genotypes.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-indels.vcf",
      " --select-type INDEL"
    ))

    ###########################################################################################
    ######### Filtering begineth 
    ###########################################################################################

    #Strict filtering
    #########################
    if (filtering.thresholds == "high") {
      # Applies filters to SNPs
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " VariantFiltration -R ", reference.path,
        " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-snps.vcf",
        " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-snps.vcf",
        " -filter \"QD<2.0\" --filter-name \"QD2\"",
        " -filter \"QUAL<30.0\" --filter-name \"QUAL30\"",
        " -filter \"SOR>3.0\" --filter-name \"SOR3\"",
        " -filter \"FS>60.0\" --filter-name \"FS60\"",
        " -filter \"MQ<40.0\" --filter-name \"MQ40\"",
        " -filter \"MQRankSum<-12.5\" --filter-name \"MQRankSum12.5\"",
        " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum8\""
      ))

      # Applies filters to indels
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " VariantFiltration -R ", reference.path,
        " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-indels.vcf",
        " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-indels.vcf",
        " -filter \"QD<2.0\" --filter-name \"QD2\"",
        " -filter \"QUAL<30.0\" --filter-name \"QUAL30\"",
        " -filter \"FS>200.0\" --filter-name \"FS200\"",
        " -filter \"ReadPosRankSum<-20.0\" --filter-name \"ReadPosRankSum20\""
      ))
    } # end high if
    
    #Medium amount of filtering
    #########################
    if (filtering.thresholds == "medium") {
      # Applies filters to SNPs
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " VariantFiltration -R ", reference.path,
        " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-snps.vcf",
        " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-snps.vcf",
        " -filter \"QD<1.0\" --filter-name \"QD2\"",
        " -filter \"QUAL<10.0\" --filter-name \"QUAL20\"",
        " -filter \"SOR>3.0\" --filter-name \"SOR3\"",
        " -filter \"FS>60.0\" --filter-name \"FS60\"",
        " -filter \"MQ<10.0\" --filter-name \"MQ40\"",
        " -filter \"MQRankSum<-12.5\" --filter-name \"MQRankSum12.5\"",
        " -filter \"ReadPosRankSum<-8.0\" --filter-name \"ReadPosRankSum8\""
      ))

      # Applies filters to indels
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " VariantFiltration -R ", reference.path,
        " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-indels.vcf",
        " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-indels.vcf",
        " -filter \"QD<2.0\" --filter-name \"QD2\"",
        " -filter \"QUAL<10.0\" --filter-name \"QUAL20\"",
        " -filter \"FS>200.0\" --filter-name \"FS200\"",
        " -filter \"ReadPosRankSum<-20.0\" --filter-name \"ReadPosRankSum20\""
      ))
    } # end high if

    #Low amount of filtering
    #########################
    if (filtering.thresholds == "low") {
      # Applies filters to SNPs
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " VariantFiltration -R ", reference.path,
        " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-snps.vcf",
        " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-snps.vcf",
      ))

      # Applies filters to indels
      system(paste0(
        gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
        " VariantFiltration -R ", reference.path,
        " -V ", output.directory, "/", sample.names[i], "/gatk4-unfiltered-indels.vcf",
        " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-indels.vcf",
      ))
    } # end high if

    ###########################################################################################
    ######### Filtering finished
    ###########################################################################################

    # Combine them into a single VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
      " SortVcf",
      " -I ", output.directory, "/", sample.names[i], "/gatk4-filtered-snps.vcf",
      " -I ", output.directory, "/", sample.names[i], "/gatk4-filtered-indels.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-filtered-genotypes.vcf"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/", sample.names[i], "/gatk4-filtered-genotypes.vcf",
      " -O ", output.directory, "/", sample.names[i], "/gatk4-final-genotypes.vcf",
      " --exclude-filtered TRUE"
    ))

    print(paste0(sample.names[i], " completed GATK4 sample genotyping!"))

  }#end i loop

  snow::stopCluster(cl)

}#end function

# END SCRIPT
