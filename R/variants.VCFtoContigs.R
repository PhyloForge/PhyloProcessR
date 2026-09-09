#' Apply filtered variants and depth rules to sample assemblies
#'
#' Produces an alternate-reference sequence from passing variants, optionally
#' encoding supported heterozygous SNP genotypes with IUPAC codes. Reference
#' bases remain outside applied variants. This is not phased haplotype output.
#'
#' @param genotype.directory Per-sample genotype directory.
#' @param mapping.directory Per-sample mapping/reference directory.
#' @param output.directory FASTA output directory.
#' @param vcf.file One of `SNP`, `Indel`, or `Both`.
#' @param consensus.sequences Produce ordinary alternate-reference output.
#' @param ambiguity.codes Produce IUPAC output.
#' @param threads Concurrent samples.
#' @param memory Total JVM heap budget in GB.
#' @param temp.directory Temporary directory.
#' @param gatk4.path GATK path/directory or NULL.
#' @param overwrite Recompute output.
#' @param quiet Suppress tool output while retaining logs.
#' @param sample.names Optional retained sample set.
#' @param depth.files Optional named depth-table paths shared between output modes.
#' @param depth.filter.mode `none`, `site`, `mean`, or `both` (default `site`).
#' @param min.site.depth Sites below this depth become N; default 1 retains 1x.
#' @param min.mean.depth Contigs below this full-span mean are removed.
#' @param max.n.proportion Optional final N-fraction cutoff.
#' @param use.base.recalibration Select recalibrated BAM when calculating depth.
#' @param samtools.path samtools path/directory or NULL.
#' @param ploidy Used only for a targeted non-diploid IUPAC warning.
#' @return Invisibly returns sample names.
#' @export
VCFtoContigs = function(genotype.directory=NULL,mapping.directory=NULL,
 output.directory="sample-contigs",vcf.file="SNP",consensus.sequences=FALSE,
 ambiguity.codes=TRUE,threads=1,memory=1,temp.directory=NULL,gatk4.path=NULL,
 overwrite=FALSE,quiet=TRUE,sample.names=NULL,depth.files=NULL,
 depth.filter.mode="site",min.site.depth=1,min.mean.depth=1,max.n.proportion=NULL,
 use.base.recalibration=FALSE,samtools.path=NULL,ploidy=2) {
  choice=match.arg(tolower(vcf.file),c("snp","indel","both")); vcf.file=c(snp="SNP",indel="Indel",both="Both")[[choice]]
  if(xor(isTRUE(consensus.sequences),isTRUE(ambiguity.codes))==FALSE) stop("Exactly one of consensus.sequences and ambiguity.codes must be TRUE.")
  settings=.validateDepthSettings(depth.filter.mode,min.site.depth,min.mean.depth,max.n.proportion)
  if(vcf.file!="SNP"&&settings$mode%in%c("site","both")) stop("Site depth masking currently requires vcf.file = 'SNP'; use mean or none for indel output.")
  if(ambiguity.codes&&ploidy!=2) warning("GATK IUPAC output does not fully encode non-diploid genotypes.")
  if(is.null(genotype.directory)||!dir.exists(genotype.directory)) stop("Genotype directory not found.")
  if(is.null(mapping.directory)||!dir.exists(mapping.directory)) stop("Mapping directory not found.")
  discovered=list.dirs(genotype.directory,recursive=FALSE,full.names=FALSE); if(is.null(sample.names)) sample.names=discovered
  if(!length(sample.names)) stop("No samples are available for FASTA conversion.")
  vname=switch(vcf.file,SNP="gatk4-final-snps.vcf",Indel="gatk4-final-indels.vcf",Both="gatk4-final-genotypes.vcf")
  vcfs=file.path(genotype.directory,sample.names,vname); refs=file.path(mapping.directory,sample.names,"index","reference.fa")
  if(any(!file.exists(vcfs))) stop("Missing selected VCF for sample(s): ",paste(sample.names[!file.exists(vcfs)],collapse=", "))
  if(any(!file.exists(refs))) stop("Missing reference for sample(s): ",paste(sample.names[!file.exists(refs)],collapse=", "))
  need.depth=settings$mode!="none"||!is.null(max.n.proportion)
  if(need.depth&&is.null(depth.files)) depth.files=calculateSampleDepth(mapping.directory,file.path(dirname(output.directory),"depth"),sample.names,use.base.recalibration,samtools.path,overwrite,quiet)
  resources=.validateResources(threads,memory,length(sample.names)); gatk=.toolCommand("gatk",gatk4.path)
  if(is.null(temp.directory)) temp.directory=tempdir(); .ensureDirectory(temp.directory)
  if(overwrite&&dir.exists(output.directory)) unlink(output.directory,recursive=TRUE); .ensureDirectory(output.directory); .ensureDirectory("logs/sample_logs")
  command=.gatkCommand(gatk,temp.directory,resources$heap.mb); mode=if(ambiguity.codes) "iupac" else "consensus"
  results=parallel::mclapply(seq_along(sample.names),function(i) tryCatch({
    s=sample.names[i]; final=file.path(output.directory,paste0(s,".fa")); marker=.stageMarker(output.directory,paste0(s,"-VCFtoContigs-",mode)); report=file.path(output.directory,paste0(s,".",mode,".depth-filter.tsv"))
    details=c(paste0("vcf=",normalizePath(vcfs[i])),paste0("mode=",settings$mode),paste0("min.site.depth=",settings$min.site.depth),paste0("min.mean.depth=",settings$min.mean.depth),paste0("max.n.proportion=",if(is.null(settings$max.n.proportion))"NULL" else settings$max.n.proportion))
    if(!overwrite&&file.exists(marker)&&file.exists(final)&&identical(readLines(marker,warn=FALSE),c("complete=true",details))) return(list(success=TRUE))
    file.remove(marker); tmp.raw=tempfile(paste0(s,"-raw-"),tmpdir=output.directory,fileext=".fa"); tmp.done=tempfile(paste0(s,"-processed-"),tmpdir=output.directory,fileext=".fa")
    log=file.path("logs/sample_logs",s,paste0("VCFtoContigs-",mode,".stderr.log")); .ensureDirectory(dirname(log))
    iupac=if(ambiguity.codes) paste("--use-iupac-sample",shQuote(s)) else ""
    .runCommand(paste(command,"FastaAlternateReferenceMaker -R",shQuote(refs[i]),"-V",shQuote(vcfs[i]),"-O",shQuote(tmp.raw),iupac),quiet,"FASTA conversion",stderr.log=log)
    if(!file.exists(tmp.raw)) stop("GATK did not create its temporary FASTA")
    seqs=Biostrings::readDNAStringSet(tmp.raw); names(seqs)=sub("^[0-9]+ ","",names(seqs))
    if(need.depth) { df=depth.files[[s]]; if(is.null(df)||!file.exists(df)) stop("Depth result missing for ",s); seqs=.filterDepthSequences(seqs,df,settings,report) }
    Biostrings::writeXStringSet(seqs,tmp.done,format="fasta",width=1000000)
    if(!file.rename(tmp.done,final)) stop("Could not publish final FASTA for ",s)
    file.remove(c(tmp.raw,paste0(tmp.raw,".fai"),sub("\\.fa$",".dict",tmp.raw)))
    writeLines(c("complete=true",details),marker); list(success=TRUE)
  },error=function(e) list(success=FALSE,message=conditionMessage(e))),mc.cores=resources$workers)
  .collectWorkers(results,sample.names,"FASTA conversion"); invisible(sample.names)
}
