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

variants.databaseSummary = function(database.directory = NULL,
                                    output.file = "database-summary",
                                    reference.path = "ref-index",
                                    threads = 1,
                                    memory = 1,
                                    resume = TRUE,
                                    overwrite = TRUE,
                                    quiet = TRUE) {

  #Debugging
  #Home directoroies
  # database.directory = "variant-calling/variant-database"
  # reference.path = "ref-index"
  # #subreference.name = "rag1"
  # output.file = "variant-calling/database-summary"
  # threads = 4 #number of threads, 8-10 is a recommended amount
  # memory = 8
  # resume = FALSE
  # quiet = FALSE
  # overwrite = TRUE

  #Quick checks
  if (is.null(database.directory) == TRUE){ stop("Please provide the directory to the variant database.") }
  if (file.exists(database.directory) == F){ stop("Database directory not found.") }

  #Creates output directory
  if (dir.exists("logs") == F){ dir.create("logs") }
  if (dir.exists(output.database) == F){ dir.create(output.database) }

  #Get loci names I guess
  loci.names = list.files(database.directory, full.names = TRUE)
  loci.names = loci.names[grep(".vcf$", loci.names)]

  if (length(loci.names) == 0){ return("no loci remain to analyze.") }

  ############################################################################################
  ########### Step 1 #########################################################################
  ##### Start up loop for each sample
  ############################################################################################

  #Sets up filtering rules to get counts
  Q20rule = TVTB::VcfFixedRules(exprs = list(qual20 = expression(QUAL >= 20) ))

  SNPrule = TVTB::VcfFixedRules(exprs = list(SNP = expression(as.character(REF) %in% c("A", "T", "G", "C") &
                                                          as.character(ALT) %in% c("A", "T", "G", "C"))))

  Indrule = TVTB::VcfFixedRules(exprs = list(INDEL = expression(Biostrings::width(REF) > Biostrings::width(ALT)) ))

  #Sets up header
  header.data = c("Locus", "TotalQual", "TotalVariants", "TotalSNPs", "TotalIndels")
  #Sets up data to collect
  collect.data = data.table::data.table(matrix(as.numeric(0), nrow = length(loci.names), ncol = length(header.data)))
  data.table::setnames(collect.data, header.data)
  collect.data[, Locus:=as.character(Locus)]

  for (i in 1:length(loci.names)){

    #Qual score for SNP missing, if needed?
    cvcf = VariantAnnotation::readVcf(file = paste0(loci.names[i]))

    #Catches empty SNPs vcf
    if (dim(cvcf)[1] == 0){
      #Collects the data
      data.table::set(collect.data, i = as.integer(i), j = match("Locus", header.data), value = loci.names[i] )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalVariants", header.data), value = 0 )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalSNPs", header.data), value = 0 )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalIndels", header.data), value = 0 )
      data.table::set(collect.data, i = as.integer(i), j = match("TotalQual", header.data), value = 0 )
      system(paste0("rm -r ", loci.names[i]))

      next
    }

    evcf = VariantAnnotation::expand(x = cvcf, row.names = TRUE)

    #Gets the counts for SNPS and indels
    q.stats = summary(TVTB::evalSeparately(Q20rule, evcf, enclos = .GlobalEnv))
    qual.20 = gsub("TRUE:", "", q.stats[2])
    qual.20 = as.numeric(gsub(" ", "", qual.20))

    all.stats = summary(TVTB::evalSeparately(SNPrule, evcf, enclos = .GlobalEnv))
    tot.snp = gsub("TRUE :", "", all.stats[3])
    tot.snp = as.numeric(gsub(" ", "", tot.snp))

    ind.stats = summary(TVTB::evalSeparately(Indrule, evcf, enclos = .GlobalEnv))
    tot.ind = gsub("TRUE :", "", ind.stats[3])
    tot.ind = as.numeric(gsub(" ", "", tot.ind))

    #Collects the data
    data.table::set(collect.data, i = as.integer(i), j = match("Locus", header.data), value = gsub(".vcf$", "", loci.names[i]) )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalQual", header.data), value = qual.20 )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalVariants", header.data), value = tot.snp + tot.ind )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalSNPs", header.data), value = tot.snp )
    data.table::set(collect.data, i = as.integer(i), j = match("TotalIndels", header.data), value = tot.ind )
    print(paste0(loci.names[i], " done!!"))

  }#end i loop

  write.table(collect.data, file = paste0(output.file, ".txt"), sep = "\t", row.names = F)

}#end functino


#########################
###### END SCRIPT
#########################
