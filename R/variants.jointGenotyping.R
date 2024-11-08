#' @title jointGenotyping
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

jointGenotyping = function(haplotype.caller.directory = "haplotype-caller",
                          output.directory = "genotype-database",
                          use.base.recalibration = FALSE,
                          save.unfiltered = TRUE,
                          save.SNPs = TRUE,
                          save.indels = TRUE,
                          save.combined = TRUE,
                          custom.SNP.QD =  NULL,
                          custom.SNP.QUAL =  NULL,
                          custom.SNP.SOR =  NULL,
                          custom.SNP.FS =  NULL,
                          custom.SNP.MQ =  NULL,
                          custom.SNP.MQRankSum =  NULL,
                          custom.SNP.ReadPosRankSum =  NULL,
                          custom.INDEL.QD =  NULL,
                          custom.INDEL.QUAL =  NULL,
                          custom.INDEL.FS =  NULL,
                          custom.INDEL.ReadPosRankSum =  NULL,
                          gatk4.path = NULL,
                          temp.directory = NULL,
                          threads = 1,
                          memory = 1,
                          overwrite = TRUE,
                          quiet = TRUE) {

 #Debugging
  # #Home directoroies
  # library(PhyloProcessR)
  # setwd("/Volumes/LaCie/Anax")

  # custom.SNP.QD <- 2
  # custom.SNP.QUAL <- 30
  # custom.SNP.SOR <- 3
  # custom.SNP.FS <- 60
  # custom.SNP.MQ <- 40
  # custom.SNP.MQRankSum <- -12.5
  # custom.SNP.ReadPosRankSum <- -8
  # custom.INDEL.QD <- 2
  # custom.INDEL.QUAL <- 30
  # custom.INDEL.FS <- 60
  # custom.INDEL.ReadPosRankSum <- -8
  # save.unfiltered = TRUE
  # save.SNPs = TRUE
  # save.indels = TRUE
  # save.combined = TRUE

  # dataset.name = "joint-genotyping"
  # haplotype.caller.directory = paste0("data-analysis/", dataset.name, "/haplotype-caller")
  # output.directory = paste0("data-analysis/", dataset.name, "/genotype-database")

  # gatk4.path <- "/Users/chutter/Bioinformatics/miniconda3/envs/PhyloProcessR/bin"
  # use.base.recalibration = FALSE
  # threads <- 4
  # memory <- 8
  # quiet <- FALSE
  # overwrite <- FALSE

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

  if (is.null(temp.directory) == TRUE){
    temp.directory = getwd()
  }


  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }

  # Sets directory and reads in  if (is.null(output.dir) == TRUE){ stop("Please provide an output directory.") }
  if (dir.exists(output.directory) == FALSE) {
    dir.create(output.directory)
    dir.create(paste0(output.directory, "/unfiltered-all"))
    dir.create(paste0(output.directory, "/unfiltered-snps"))
    dir.create(paste0(output.directory, "/unfiltered-indels"))
    dir.create(paste0(output.directory, "/filtered-all"))
    dir.create(paste0(output.directory, "/filtered-snps"))
    dir.create(paste0(output.directory, "/filtered-indels"))

  } else {
    if (overwrite == TRUE) {
      system(paste0("rm -r ", output.directory))
      dir.create(output.directory)
      dir.create(paste0(output.directory, "/unfiltered-all"))
      dir.create(paste0(output.directory, "/unfiltered-snps"))
      dir.create(paste0(output.directory, "/unfiltered-indels"))
      dir.create(paste0(output.directory, "/filtered-all"))
      dir.create(paste0(output.directory, "/filtered-snps"))
      dir.create(paste0(output.directory, "/filtered-indels"))
    }
  } # end else

  SNP.QD.string <- paste0(" -filter \"QD<", format(custom.SNP.QD, nsmall = 1), "\" --filter-name \"QD\"")
  SNP.QUAL.string <- paste0(" -filter \"QUAL<", format(custom.SNP.QUAL, nsmall = 1), "\" --filter-name \"QUAL\"")
  SNP.SOR.string <- paste0(" -filter \"SOR>", format(custom.SNP.SOR, nsmall = 1), "\" --filter-name \"SOR\"")
  SNP.FS.string <- paste0(" -filter \"FS>", format(custom.SNP.FS, nsmall = 1), "\" --filter-name \"FS\"")
  SNP.MQ.string <- paste0(" -filter \"MQ<", format(custom.SNP.MQ, nsmall = 1), "\" --filter-name \"MQ\"")
  SNP.MQRankSum.string <- paste0(" -filter \"MQRankSum<", format(custom.SNP.MQRankSum, nsmall = 1), "\" --filter-name \"MQRankSum\"")
  SNP.ReadPosRankSum.string <- paste0(" -filter \"ReadPosRankSum<", format(custom.SNP.ReadPosRankSum, nsmall = 1), "\" --filter-name \"ReadPosRankSum\"")

  IN.QD.string <- paste0(" -filter \"QD<", format(custom.INDEL.QD, nsmall = 1), "\" --filter-name \"QD\"")
  IN.QUAL.string <- paste0(" -filter \"QUAL<", format(custom.INDEL.QUAL, nsmall = 1), "\" --filter-name \"QUAL\"")
  IN.FS.string <- paste0(" -filter \"FS>", format(custom.INDEL.FS, nsmall = 1), "\" --filter-name \"FS\"")
  IN.ReadPosRankSum.string <- paste0(" -filter \"ReadPosRankSum<", format(custom.INDEL.ReadPosRankSum, nsmall = 1), "\" --filter-name \"ReadPosRankSum\"")

  # Get loci names I guess
  reference.path = "index/reference.fa"
  reference.seq = Biostrings::readDNAStringSet((paste0("index/reference.fa")))
  loci.names = names(reference.seq)

  # Resumes file download
  if (overwrite == FALSE) {
    done.files <- list.files(paste0(output.directory, "/filtered-all"), full.names = T, recursive = T)
    done.names <- gsub(".vcf$", "", done.files)
    done.names <- gsub(".*\\/", "", done.names)
    loci.names <- loci.names[!loci.names %in% done.names]
  }

  # Get multifile databases together
  sample.names <- list.files(haplotype.caller.directory, recursive = TRUE, full.names = TRUE)

  if (length(sample.names) == 0) {
    return("no samples available to analyze.")
  }

  if (length(sample.names) == 0) {
    return("no samples available to analyze.")
  }

  if (use.base.recalibration == TRUE) {
    sample.names = sample.names[grep("gatk4-bqsr-haplotype-caller.g.vcf.gz$", sample.names)]
    if (length(sample.names) == 0) {
      stop("the haplotype caller file for use.base.recalibration does not exist.")
    }
  }

  if (use.base.recalibration == FALSE) {
    sample.names = sample.names[grep("gatk4-haplotype-caller.g.vcf.gz$", sample.names)]
    if (length(sample.names) == 0) {
      stop("the haplotype caller file does not exist.")
    }
  }

  #Gathers input file names together
  sample.datasets = paste0("-V ", sample.names)
  vcf.files = paste0(sample.datasets, collapse = " ")

  # Sets up multiprocessing
  cl = snow::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  mem.cl = floor(memory / threads)

  # Loops through each locus and does operations on them
  foreach::foreach(i = seq_along(loci.names), .packages = c("foreach")) %dopar% {
    # Loops through each locus and does operations on them
    # for (i in 1:length(loci.names)){

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " GenomicsDBImport ", vcf.files,
      " --genomicsdb-workspace-path ", output.directory, "/", loci.names[i],
      " --intervals ", loci.names[i]
    ))

    # Genotype haplotype caller results
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " GenotypeGVCFs -R ", reference.path,
      " -V gendb://", output.directory, "/", loci.names[i],
      " --use-new-qual-calculator true",
      " -O ", output.directory, "/unfiltered-all/", loci.names[i], ".vcf"
    ))

    # Selects only the SNPs from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-all/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-snps/", loci.names[i], ".vcf",
      " --select-type SNP"
    ))

    # Selects only the indels from the VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-all/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-indels/", loci.names[i], ".vcf",
      " --select-type INDEL"
    ))

    ###########################################################################################
    ######### Filtering begineth
    ###########################################################################################

    # Custom filtering
    #########################
    # Applies filters to SNPs
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", output.directory, "/unfiltered-snps/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-snps/", loci.names[i], "_filter.vcf",
      SNP.QD.string,
      SNP.QUAL.string,
      SNP.SOR.string,
      SNP.FS.string,
      SNP.MQ.string,
      SNP.MQRankSum.string,
      SNP.ReadPosRankSum.string
    ))

    # Applies filters to indels
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " VariantFiltration -R ", reference.path,
      " -V ", output.directory, "/unfiltered-indels/", loci.names[i], ".vcf",
      " -O ", output.directory, "/unfiltered-indels/", loci.names[i], "_filter.vcf",
      IN.QD.string,
      IN.QUAL.string,
      IN.FS.string,
      IN.ReadPosRankSum.string
    ))

    ###########################################################################################
    ######### Filtering finished
    ###########################################################################################

    # Combine them into a single VCF
    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SortVcf",
      " -I ", output.directory, "/unfiltered-snps/", loci.names[i], "_filter.vcf",
      " -I ", output.directory, "/unfiltered-indels/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/unfiltered-all/", loci.names[i], "_filter.vcf"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-all/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/filtered-all/", loci.names[i], ".vcf",
      " --exclude-filtered TRUE"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-snps/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/filtered-snps/", loci.names[i], ".vcf",
      " --exclude-filtered TRUE"
    ))

    system(paste0(
      gatk4.path, "gatk --java-options \"-Djava.io.tmpdir=", temp.directory, " -Xmx", memory, "G\"",
      " SelectVariants",
      " -V ", output.directory, "/unfiltered-indels/", loci.names[i], "_filter.vcf",
      " -O ", output.directory, "/filtered-indels/", loci.names[i], ".vcf",
      " --exclude-filtered TRUE"
    ))

    print(paste0(loci.names[i], " completed GATK4 joint sample genotyping!"))
    system(paste0("rm -r ", output.directory, "/", loci.names[i]))

  }#end i loop

  snow::stopCluster(cl)

  #Datasets to save
  #Save SNPs
  if (save.SNPs == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered-snps ", output.directory, "/filtered-snps"))
  }

  #Save indels
  if (save.indels == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered-indels ", output.directory, "/filtered-indels"))
  }

  #Save all
  if (save.combined == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered-all ", output.directory, "/filtered-all"))
  }

  #Save unfiltered
  if (save.unfiltered == FALSE) {
    system(paste0("rm -r ", output.directory, "/unfiltered*"))
  }

}#end function

# END SCRIPT
