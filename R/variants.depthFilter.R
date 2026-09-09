# Internal depth helpers for workflow 3.
.validateDepthSettings = function(depth.filter.mode="site", min.site.depth=1,
                                  min.mean.depth=1, max.n.proportion=NULL) {
  mode=match.arg(tolower(depth.filter.mode),c("none","site","mean","both"))
  for(n in c("min.site.depth","min.mean.depth")) { x=get(n); if(length(x)!=1||!is.finite(x)||x<0) stop(n," must be one non-negative finite value.") }
  if(!is.null(max.n.proportion)&&(length(max.n.proportion)!=1||!is.finite(max.n.proportion)||max.n.proportion<0||max.n.proportion>1)) stop("max.n.proportion must be NULL or one value from 0 to 1.")
  list(mode=mode,min.site.depth=min.site.depth,min.mean.depth=min.mean.depth,max.n.proportion=max.n.proportion)
}

#' Calculate full-span per-base depth for workflow 3
#' @param mapping.directory Mapped sample directory.
#' @param output.directory Owned depth-table directory.
#' @param sample.names Samples to scan.
#' @param use.base.recalibration Use recalibrated BAMs.
#' @param samtools.path samtools path/directory or NULL.
#' @param overwrite Recompute cached depth.
#' @param quiet Suppress samtools output while retaining logs.
#' @return Named paths to depth tables.
#' @export
calculateSampleDepth = function(mapping.directory,output.directory,sample.names=NULL,
 use.base.recalibration=FALSE,samtools.path=NULL,overwrite=FALSE,quiet=TRUE) {
  if(!dir.exists(mapping.directory)) stop("Mapping directory not found.")
  if(is.null(sample.names)) sample.names=list.dirs(mapping.directory,recursive=FALSE,full.names=FALSE)
  if(!length(sample.names)) stop("No samples are available for depth calculation.")
  samtools=.toolCommand("samtools",samtools.path); .ensureDirectory(output.directory); .ensureDirectory("logs/sample_logs")
  paths=setNames(file.path(output.directory,paste0(sample.names,".depth.tsv")),sample.names)
  for(s in sample.names) {
    bam=.selectedSampleBam(mapping.directory,s,use.base.recalibration); ref=file.path(mapping.directory,s,"index","reference.fa"); fai=paste0(ref,".fai")
    if(!file.exists(fai)) stop("Reference index not found for ",s)
    record=paste0(paths[s],".complete"); info=file.info(bam)
    signature=c(paste0("bam=",normalizePath(bam)),paste0("size=",info$size),paste0("mtime=",as.numeric(info$mtime)),paste0("reference=",normalizePath(ref)),"flags=-aa -s -q 0 -Q 0 -G 0x800")
    reusable=!overwrite&&file.exists(paths[s])&&file.info(paths[s])$size>0&&file.exists(record)&&identical(readLines(record,warn=FALSE),signature)
    if(reusable) next
    file.remove(c(paths[s],record)); log=file.path("logs/sample_logs",s,"depth.stderr.log"); .ensureDirectory(dirname(log))
    cmd=paste(samtools,"depth -aa -s -q 0 -Q 0 -G 0x800",shQuote(bam),">",shQuote(paths[s]))
    .runCommand(cmd,quiet,"samtools depth",keep.stdout=TRUE,stderr.log=log)
    if(!file.exists(paths[s])||file.info(paths[s])$size==0) stop("Depth output is empty for ",s)
    # Validate shape, contigs, full lengths, and numeric coordinates without loading the BAM.
    d=data.table::fread(paths[s],header=FALSE,select=1:3,col.names=c("contig","position","depth"),showProgress=FALSE)
    ix=data.table::fread(fai,header=FALSE,select=1:2,col.names=c("contig","length"),showProgress=FALSE)
    if(any(!d$contig%in%ix$contig)||any(!is.finite(d$position))||any(!is.finite(d$depth))||any(d$depth<0)) stop("Invalid or truncated depth table for ",s)
    n.by = table(d$contig)
    counts = data.frame(contig = names(n.by), n = as.integer(n.by),
                        min = as.numeric(tapply(d$position, d$contig, min)),
                        max = as.numeric(tapply(d$position, d$contig, max)))
    chk=merge(ix,counts,by="contig",all.x=TRUE)
    if(any(is.na(chk$n))||any(chk$n!=chk$length)||any(chk$min!=1)||any(chk$max!=chk$length)) stop("Depth table does not cover every indexed reference position for ",s)
    writeLines(signature,record)
  }
  paths
}

.filterDepthSequences = function(seqs,depth.file,settings,report.file) {
  d=data.table::fread(depth.file,header=FALSE,col.names=c("contig","position","depth"),showProgress=FALSE)
  original.names=names(seqs); normalized=sub("^[0-9]+ ","",original.names)
  if(anyDuplicated(normalized)||!setequal(normalized,unique(d$contig))) stop("FASTA and depth contig names do not agree.")
  names(seqs)=normalized; rows=vector("list",length(seqs)); keep=rep(TRUE,length(seqs))
  for(i in seq_along(seqs)) {
    x=d[d$contig==names(seqs)[i]]; len=Biostrings::width(seqs[i]); if(nrow(x)!=len) stop("Depth length mismatch for contig ",names(seqs)[i])
    old.n=as.integer(Biostrings::letterFrequency(seqs[i],"N")); mean.depth=mean(x$depth); below=sum(x$depth<settings$min.site.depth); newly=0L
    if(settings$mode%in%c("site","both")) { pos=x$position[x$depth<settings$min.site.depth]; before=as.character(seqs[i]); if(length(pos)) seqs[i]=Biostrings::replaceAt(seqs[i],IRanges::IRanges(pos,width=1),Biostrings::DNAStringSet(rep("N",length(pos)))); newly=sum(substring(before,pos,pos)!="N") }
    final.n=as.integer(Biostrings::letterFrequency(seqs[i],"N")); prop=if(len) final.n/len else 1; reasons=character()
    if(settings$mode%in%c("mean","both")&&mean.depth<settings$min.mean.depth) reasons=c(reasons,"mean_depth")
    filtering=settings$mode!="none"||!is.null(settings$max.n.proportion)
    if(filtering&&(len==0||final.n==len)) reasons=c(reasons,"no_sequence")
    if(!is.null(settings$max.n.proportion)&&prop>settings$max.n.proportion) reasons=c(reasons,"excess_N")
    keep[i]=!length(reasons); rows[[i]]=data.frame(contig=names(seqs)[i],length=len,mean_depth=mean.depth,bases_below_threshold=below,preexisting_N=old.n,newly_masked=newly,final_N_proportion=prop,retained=keep[i],reason=paste(unique(reasons),collapse=";"))
  }
  data.table::fwrite(data.table::rbindlist(rows),report.file,sep="\t")
  seqs[keep]
}
