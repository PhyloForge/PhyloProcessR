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

variants.databaseSummary = function(bam.directory = NULL,
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
  loci.names = list.files(db.dir, full.names = TRUE)
  loci.names = loci.names[grep(".vcf$", loci.names)]

  #Get multifile databases together
  sample.names = list.dirs(bam.directory, recursive = F, full.names = F)
  all.files = list.files(bam.directory, recursive = T, full.names = T)
  vcf.files = all.files[grep("gatk4-haplotype-caller.g.vcf.gz$", all.files)]

  if (length(sample.names) == 0){ return("no samples remain to analyze.") }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  #Sets up filtering rules to get counts
  SNPrule = VcfFixedRules(exprs = list(qual20 = expression(QUAL >= 20),
                                       SNP = expression(as.character(REF) %in% c("A", "T", "G", "C") &
                                                          as.character(ALT) %in% c("A", "T", "G", "C"))))

  Indrule = VcfFixedRules(exprs = list(INDEL = expression(Biostrings::width(REF) > Biostrings::width(ALT)) ))

  #Sets up header
  header.data = c("Locus", "TotalVariants", "TotalSNPs", "TotalIndels", "TotalQual", "TotalPass")
  #Sets up data to collect
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(loci.names), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, Locus:=as.character(Locus)]

  for (i in 1:length(loci.names)){

    #Qual score for SNP missing, if needed?
    cvcf = VariantAnnotation::readVcf(file = paste0(loci.names[i], ".vcf"))

    #Catches empty SNPs vcf
    if (dim(cvcf)[1] == 0){
      #Collects the data
      data.table::set(collect.data, i = as.integer(i), j = match("Locus", header.data), value = loci.names[i] )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalVariants", header.data), value = 0 )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalSNPs", header.data), value = 0 )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalIndels", header.data), value = 0 )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalQual", header.data), value = 0 )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalPass", header.data), value = 0 )
      system(paste0("rm -r ", loci.names[i]))

      next
    }

    evcf = VariantAnnotation::expand(x = cvcf, row.names = TRUE)

    #Gets the counts for SNPS and indels
    all.stats = summary(VariantAnnotation::evalSeparately(SNPrule, evcf, enclos = .GlobalEnv))
    ind.stats = summary(VariantAnnotation::evalSeparately(Indrule, evcf, enclos = .GlobalEnv))

    #Collects the data
    data.table::set(collect.data, i = as.integer(i), j = match("Locus", header.data), value = loci.names[i] )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalVariants", header.data), value = as.numeric(all.stats[1]) )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalSNPs", header.data), value = as.numeric(all.stats[3]) )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalIndels", header.data), value = as.numeric(ind.stats[2]) )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalQual", header.data), value = as.numeric(all.stats[2]) )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalPass", header.data), value = as.numeric(all.stats[4]) )

    print(paste0(loci.names[i], " done!!"))

  }#end i loop

  write.table(collect.data, file = "Locus_SNP_Summary.txt", sep = "\t", row.names = F)

}#end functino


#########################
###### END SCRIPT
#########################
