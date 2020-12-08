#Libraries
library(rdrop2)


#################################################
## Step 0: Data download. Optional step.
##################

### Example usage
sample.spread = "/home/c111h652/scratch/MitoGenomes/Mitogenome_study.csv"
out.dir = "/home/c111h652/scratch/MitoGenomes/raw-reads-frogs"
dropbox.dir = "/Research/3_Sequence-Database/Raw-Reads"
dropbox.tok = "/home/c111h652/dropbox-token.RDS"

#Run download function
dropboxDownload(sample.spreadsheet = sample.spread,
                dropbox.directory = dropbox.dir,
                dropbox.token = dropbox.tok,
                out.directory = out.dir,
                overwrite = TRUE)


#################################################
## Step 1: Preprocess reads with read processing function
##################






#################################################
## Step 2: Clean out contamination
##################






#################################################
## Step 3: iteratively assemble mt-genomes
##################







#################################################
## Step 4: do all the alignment stuff
##################






###############################################################################
###############################################################################
######################           FUNCTIONS            #########################
###############################################################################
###############################################################################

#Consensus function
make.consensus =  function (input.alignment, method = c("majority", "threshold", "IUPAC",
                                                        "profile"), threshold = 0.6, warn.non.IUPAC = FALSE, type = c("DNA", "RNA")) {

  #input.alignment<-trimmed
  #Converts alignment to matrix of characters to be used
  new.align<-strsplit(as.character(input.alignment), "")
  align.in<-matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)

  #Does based on method
  method <- match.arg(method)

  if (method == "IUPAC") {
    type <- match.arg(type)
    res <- apply(align.in, 2, bma, warn.non.IUPAC = warn.non.IUPAC,
                 type = type)
    names(res) <- NULL
  }
  if (method == "majority") {
    majority <- function(x) names(which.max(table(x)))
    res <- apply(align.in, 2, majority)
    names(res) <- NULL
  }
  if (method == "profile") {
    obsvalue <- levels(factor(align.in))
    nrow <- length(obsvalue)
    row.names(align.in) <- NULL
    res <- apply(matali, 2, function(x) table(factor(x, levels = obsvalue)))
  }
  if (method == "threshold") {
    profile <- consensus(align.in, method = "profile")
    profile.rf <- apply(profile, 2, function(x) x/sum(x))
    res <- rownames(profile.rf)[apply(profile.rf, 2, which.max)]
    res <- ifelse(apply(profile.rf, 2, max) >= threshold,
                  res, NA)
    names(res) <- NULL
  }

  out.consensus<-DNAStringSet(paste0(res, collapse = ""))
  names(out.consensus)<-"Consensus_Sequence"
  return(out.consensus)
}
#Slice function
slice.trim = function(input.align, slice.size.bp = 100, threshold = 0.45){

  #makes consensus sequence for comparison
  #input.align<-align
  input.con<-make.consensus(input.align, method = "majority")
  names(input.con)<-"Reference_Locus"

  comb.align<-append(input.align, input.con)

  #Gets slice information ready
  slice.no<-ceiling(max(width(input.align))/slice.size.bp)
  slice.start<-1
  slice.end<-slice.size.bp

  #checks to see if its out of bounds
  if (slice.end > max(width(input.align))){
    slice.end<-max(width(input.align))
  }#end if check
  output.align<-DNAStringSet()
  for (x in 1:slice.no){

    #Slice alignment into number of slices
    sliced.align<-subseq(comb.align, start = slice.start, end = slice.end)
    #Checks for badly aligned sequences
    bad.align<-pairwise.inf.sites(sliced.align, "Reference_Locus")
    #Remove bad sequence chunks
    rem.seqs<-bad.align[bad.align >= threshold]
    good.align<-sliced.align[!names(sliced.align) %in% names(rem.seqs)]
    #Makes replacement gap seqs for the bad ones
    blank.align<-DNAStringSet()
    if (length(rem.seqs) != 0){
      for (y in 1:length(rem.seqs)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", slice.end-slice.start+1), collapse = "")) )
      }
      names(blank.align)<-names(rem.seqs)
    }#end rem seqs if

    #Saves the slices and cats
    save.slice<-append(good.align, blank.align)
    save.slice<-save.slice[order(names(save.slice))]
    save.names<-names(save.slice)
    output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
    names(output.align)<-save.names

    #Gets new start and stop
    slice.start<-slice.start+100
    slice.end<-slice.end+100
    #checks to see if its out of bounds
    if (slice.end > max(width(input.align))){
      slice.end<-max(width(input.align))
      if (slice.end-slice.start <= 25){ break } else {
        save.slice<-subseq(comb.align, start = slice.start, end = slice.end)
        save.slice<-save.slice[order(names(save.slice))]
        save.names<-names(save.slice)
        output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
        names(output.align)<-save.names
        break
      }
    }#end if
  }#end x loop

  #Removes reference
  output.align<-output.align[names(output.align) != "Reference_Locus"]
  #removes gap only taxa
  str.splitted<-strsplit(as.character(output.align), "")
  x.align<-as.matrix(as.DNAbin(str.splitted) )
  len.temp<-as.character(as.list(x.align))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= 20]
  return.align<-output.align[!names(output.align) %in% names(spp.rem)]
  return(return.align)
}#end FUNCTION

pairwise.inf.sites2 = function(x, y) {
  #Sets up data, puts ref seq first
  new.align<-x
  new.align[new.align == "n"]<-"-"
  new.align[is.na(new.align) == T]<-"-"
  ref<-new.align[rownames(new.align) == y,]
  summary.data<-c()
  all.pars<-c()
  all.over<-c()
  for (z in 1:nrow(new.align)) {
    #Site counter
    pars<-0
    overlap<-0
    tar<-new.align[z,]
    combined<-matrix(NA_character_, ncol = max(length(ref), length(tar)), nrow =2)
    combined[1,]<-ref
    combined[2,]<-tar
    for (k in 1:ncol(combined)) {
      #Pulls out column of data
      seq.col<-vector("character", length = nrow(combined))
      seq.col<-combined[,k]
      #not equal to -
      f.char<-seq.col[seq.col != '-']
      #don't count missing seq
      if (length(f.char) <= 1) { next }

      if (length(f.char) >= 2){
        overlap<-overlap+1
        if (f.char[1] != f.char [2]) { pars<-pars+1 }
      }#end if
    }#ends informative sites loop
    all.pars<-append(all.pars, pars)
    all.over<-append(all.over, overlap)
  }# ends seq loop
  #Summarizes and returns data
  summary.data<-all.pars/all.over
  summary.data[is.nan(summary.data)]<-0
  names(summary.data)<-rownames(new.align)
  return(summary.data)
}

pairwise.inf.sites = function(x, y) {
  #Alignment should be DNAStringSet
  # x<-m.align
  # y<-"Reference_Locus"
  temp.align<-strsplit(as.character(x), "")
  mat.align<-lapply(temp.align, tolower)
  m.align<-as.matrix(as.DNAbin(mat.align))

  #Filters out weirdly divergent sequences
  new.align<-as.character(m.align)
  new.align[new.align == "n"]<-"-"
  new.align[is.na(new.align) == T]<-"-"
  ref<-new.align[rownames(new.align) == y,]
  summary.data<-c()
  all.pars<-c()
  all.over<-c()
  for (z in 1:nrow(new.align)) {
    #Site counter
    pars<-0
    overlap<-0
    tar<-new.align[z,]
    combined<-matrix(NA_character_, ncol = max(length(ref), length(tar)), nrow =2)
    combined[1,]<-ref
    combined[2,]<-tar
    for (k in 1:ncol(combined)) {
      #Pulls out column of data
      seq.col<-vector("character", length = nrow(combined))
      seq.col<-combined[,k]
      #not equal to -
      f.char<-seq.col[seq.col != '-']
      #don't count missing seq
      if (length(f.char) <= 1) { next }

      if (length(f.char) >= 2){
        overlap<-overlap+1
        if (f.char[1] != f.char [2]) { pars<-pars+1 }
      }#end if
    }#ends informative sites loop
    all.pars<-append(all.pars, pars)
    all.over<-append(all.over, overlap)
  }# ends seq loop
  #Summarizes and returns data
  summary.data<-all.pars/all.over
  summary.data[is.nan(summary.data)]<-0
  names(summary.data)<-rownames(new.align)
  return(summary.data)
}

write.phy = function (x, file = "", interleave = FALSE, strict = FALSE){
  str2cha <- function(x) {
    unlist(strsplit(x, ""))
  }
  datatype <- ifelse(is.numeric(x[1, 1]), "continuous", "nc")
  ntax <- nrow(x)
  nchar <- ncol(x)
  taxnames <- rownames(x)
  if (strict) {
    taxnames <- substring(taxnames, 1, truncate)
    missing <- 10 - unlist(lapply(strsplit(taxnames, ""),
                                  length))
    for (i in seq(along = taxnames)) taxnames[i] <- paste(taxnames[i],
                                                          paste(rep("*", missing[i]), collapse = ""), sep = "")
    if (any(duplicated(taxnames)))
      cat("WARNING: Truncation of taxon names created",
          "identical strings.")
  }
  else {
    xx <- nchar(taxnames)
    diff <- max(xx) - xx + 3
    for (i in 1:ntax) taxnames[i] <- paste(taxnames[i], paste(rep(" ",
                                                                  diff[i]), collapse = ""), sep = "")
  }
  if (!interleave)
    interleave <- nchar
  nbpart <- ceiling(nchar/interleave)
  pt <- matrix(nrow = nbpart, ncol = 2)
  pt[1, ] <- c(1, interleave)
  if (nbpart > 1)
    for (i in 2:(dim(pt)[1])) {
      pt[i, ] <- c(pt[i - 1, 2] + 1, pt[i - 1, 2] + interleave)
      pt[nbpart, 2] <- nchar
    }
  phy <- paste(ntax, nchar)
  for (i in seq(along = pt[, 1])) {
    sm <- as.character(x[, pt[i, 1]:pt[i, 2]])
    if (is.null(dim(sm)))
      sm <- as.matrix(sm, ncol = 1)
    sm <- apply(sm, 1, paste, collapse = "")
    if (i == 1)
      sm <- paste(taxnames, sm)
    if (i < max(seq(along = pt[, 1])))
      sm <- c(sm, "")
    phy <- c(phy, sm)
  }
  if (file == "") {
    cat(phy, sep = "\n")
  }
  else {
    write(phy, file = file)
  }
}

run.mafft = function(unaligned.contigs, add.contigs = NULL, algorithm = "localpair", rev.dir = T, save.name = NULL, threads = 6, delete.files = T){

  #save.name<-locus.save.name
  #algorithm = "localpair"
  # unaligned.contigs<-intron.align
  save.contigs<-as.list(as.character(unaligned.contigs))
  if (is.null(save.name) == T) { save.name<-paste(sample(LETTERS, 5, replace = T), collapse = "")}
  if (rev.dir == T){ adjust.direction<-"--adjustdirection" } else { adjust.direction<-"" }

  #Adds a sequence into the alignment. Saves much computation.
  if (algorithm == "add"){
    #Saves to folder to run with mafft
    write.fasta(sequences = save.contigs, names = names(save.contigs),
                paste0(save.name, ".fa"), nbchar = 1000000, as.string = T)

    #Saves to folder to run with mafft
    add.save<-as.list(as.character(add.contigs))
    write.fasta(sequences = add.save, names = names(add.save),
                "add_sequences.fa", nbchar = 1000000, as.string = T)

    #Runs MAFFT to align
    system(paste0(mafft.path, " --",algorithm, " add_sequences.fa --maxiterate 1000 ", adjust.direction,
                  " ", save.name, ".fa > ", save.name, "_align.fa"), ignore.stderr = T)

    alignment<-scanFa(FaFile(paste(save.name, "_align.fa", sep = "")))   # loads up fasta file
    unlink(paste(save.name, ".fa", sep = ""))
    unlink("add_sequences.fa")

  }#end -add

  #Does Regular MAFFT Local Pair
  if (algorithm == "localpair"){
    #Saves to folder to run with mafft
    write.fasta(sequences = save.contigs, names = names(save.contigs),
                paste(save.name, ".fa", sep = ""), nbchar = 1000000, as.string = T)

    #Runs MAFFT to align
    system(paste0(mafft.path, " --",algorithm, " --maxiterate 1000 ", adjust.direction, " --quiet --op 3 --ep 0.123",
                  " --thread ", threads, " ", save.name, ".fa > ", save.name, "_align.fa"), ignore.stderr = T)

    alignment<-scanFa(FaFile(paste(save.name, "_align.fa", sep = "")))   # loads up fasta file
    unlink(paste(save.name, ".fa", sep = ""))
  }#end local pair

  if (delete.files == T){
    unlink(paste(save.name, "_align.fa", sep = ""))
    return(alignment)
  } else { return(alignment) }
}#function end


find.codon = function(input.seq, genetic.code = "nuclear", reverse = T){
  #Sets up data
  #reverse = T
  #genetic.code = "mitochondrial"
  #input.seq = trimmed[j]

  #Can add custom codes in ehre
  if (genetic.code == "mitochondrial"){

    vert.mit = c("TAA", "TAG", "AGA", "AGG")
    gen.code = paste0(vert.mit, collapse = "|")

  }#end mito

  if (genetic.code == "nuclear"){

    vert.nuc = c("TAA", "TGA", "TAG")
    gen.code = paste0(vert.nuc, collapse = "|")

  }#end nuc

  if (genetic.code == "custom"){

    vert.nuc = c("TAA", "TGA", "TAG")
    gen.code = paste0(vert.nuc, collapse = "|")

  }#end nuc

  #Gets codon stuff
  for.seq = as.character(input.seq)
  codons = data.frame(str_locate_all(for.seq, gen.code)[[1]])
  codon.table = c()
  #Forward frame 1
  temp.codons = codons[(codons$start+2) %% 3 == 0,]
  if (nrow(temp.codons) != 0){
    codon.table = rbind(codon.table, data.frame(frame = "F1", temp.codons))
  } else { codon.table = rbind(codon.table, data.frame(frame = "F1", start = 0, end = 0)) }

  #Forward frame 2
  temp.codons = codons[(codons$start+1) %% 3 == 0,]
  if (nrow(temp.codons) != 0){
    codon.table = rbind(codon.table, data.frame(frame = "F2", temp.codons))
  } else { codon.table = rbind(codon.table, data.frame(frame = "F2", start = 0, end = 0)) }

  #Forward frame 3
  temp.codons = codons[(codons$start) %% 3 == 0,]
  if (nrow(temp.codons) != 0){
    codon.table = rbind(codon.table, data.frame(frame = "F3", temp.codons))
  } else { codon.table = rbind(codon.table, data.frame(frame = "F3", start = 0, end = 0)) }


  if (reverse == T){
    #Reverses seq
    rev.seq = as.character(reverseComplement(input.seq))
    #Gets codon stuff
    codons = data.frame(str_locate_all(rev.seq, gen.code)[[1]])
    #Forward frame 1
    temp.codons = codons[(codons$start+2) %% 3 == 0,]
    if (nrow(temp.codons) != 0){
      codon.table = rbind(codon.table, data.frame(frame = "R1", temp.codons))
    } else { codon.table = rbind(codon.table, data.frame(frame = "R1", start = 0, end = 0)) }

    #Forward frame 2
    temp.codons = codons[(codons$start+1) %% 3 == 0,]
    if (nrow(temp.codons) != 0){
      codon.table = rbind(codon.table, data.frame(frame = "R2", temp.codons))
    } else { codon.table = rbind(codon.table, data.frame(frame = "R2", start = 0, end = 0)) }

    #Forward frame 3
    temp.codons = codons[(codons$start) %% 3 == 0,]
    if (nrow(temp.codons) != 0){
      codon.table = rbind(codon.table, data.frame(frame = "R3", temp.codons))
    } else { codon.table = rbind(codon.table, data.frame(frame = "R3", start = 0, end = 0)) }
  }#end reverse i

  return(codon.table)

}#end function

find.orf = function(input.seq, min.size = 50, mitochondrial = F){

  #Sets up data
  #mitochondrial = T
  #min.size = 50
  #input.seq = trimmed[j]

  c.table = find.codon(input.seq, genetic.code = "mitochondrial", reverse = T)

  frames = unique(c.table$frame)
  orf.frame = data.frame()
  for (x in 1:length(frames)){
    temp.codon = c.table[c.table$frame == frames[x],]
    temp.codon = temp.codon[order(temp.codon$start),]

    if (temp.codon$start[1] == 0){
      temp.start = as.numeric(gsub("F|R", "", temp.codon$frame))
      add.frame = data.frame(FrameStart = temp.start, FrameEnd = width(input.seq),
                             Size = (width(input.seq)-temp.start)+1, Frame = frames[x])
      orf.frame = rbind(orf.frame, add.frame)
      next
    }
    #Goes through each of the given directions codons and converts to frame ranges
    temp.frame = data.frame()
    for (y in 1:(nrow(temp.codon)+1)){
      #First y the start is 1, otherwise take from previous end
      if (y == 1){ frame.start = as.numeric(gsub("F|R", "", temp.codon$frame[y])) } else { frame.start<-temp.frame$FrameEnd[y-1]+4 }

      #Gets end by subtracting from the codon start
      frame.end = temp.codon$start[y]-1
      temp.frame = rbind(temp.frame, data.frame(FrameStart = frame.start, FrameEnd = frame.end))
    } # end y loop

    temp.frame$FrameEnd[nrow(temp.frame)]<-width(input.seq)

    #Adds all the data together
    add.frame = cbind(temp.frame, Size = (temp.frame$FrameEnd-temp.frame$FrameStart)+1, Frame = frames[x])
    orf.frame = rbind(orf.frame, add.frame)

  } #end x loop

  orf.frame = orf.frame[orf.frame$Size >= min.size,]
  return(orf.frame)

}# END FUNCTION

trim.ends<-function (x, min.n.seq = 4, codon.trim = T){
  #Converts DNAStringSet to something usable
  #x<-trimmed
  #min.n.seq<-4
  new.align<-strsplit(as.character(x), "")
  mat.align<-lapply(new.align, tolower)
  x<-as.matrix(as.DNAbin(mat.align))

  if (!inherits(x, "DNAbin")) {
    stop("'x' is not of class 'DNAbin'")
  }
  if (!is.matrix(x)) {
    stop("'x' must be a matrix")
  }

  replaceWithN <- function(x) {
    id <- x == as.raw(4)
    if (length(id) > 0 & any(id[c(1, length(id))])) {
      id <- which(id)
      getIndex <- function(x) {
        for (i in seq_along(id) - 1) {
          if (any(id[1:(i + 1)] != (1:(i + 1))))
            break
        }
        id <- rev(id)
        jj <- head(id, 1)
        j <- jj - 1
        for (k in seq_along(id)[-1]) {
          if (any(id[1:k] != (jj:j)))
            break
          j <- j - 1
        }
        j <- j + 1
        id <- c(0:i, j:jj)
        id[id != 0]
      }
      id <- getIndex(id)
      x[id] <- as.raw(240)
    }
    return(x)
  }
  #Does stuff
  x <- t(apply(x, 1, replaceWithN))
  class(x) <- "DNAbin"
  b <- as.raw(c(136, 40, 72, 24))
  percentInformation <- function(x, b) {
    length(which(x %in% b))
  }
  m <- apply(x, 2, percentInformation, b)
  if (max(m) < min.n.seq) {
    return(0)
  }#end if
  m <- range(which(m >= min.n.seq))

  #Forward Frame 2
  if (codon.trim == T){
    if ((m[1]-1) %% 3 == 0){ m[1]<-m[1] }
    if ((m[1]-1) %% 3 == 1){ m[1]<-m[1]+2 }
    if ((m[1]-1) %% 3 == 2){ m[1]<-m[1]+1 }
  }

  m <- seq(from = m[1], to = m[2])
  x2 <- as.matrix(x[, m])
  #Converts back
  save.names<-rownames(x2)

  #Removes N end gaps
  x3<-as.list(data.frame(t(as.character(x2))))
  for (y in 1:length(x3)){
    #Starts from the beginning and end to fill in end gaps
    for (q in 1:length(x3[[y]])){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
    for (q in length(x3[[y]]):1){ if (x3[[y]][q] == "n"){ x3[[y]][q]<-"-" } else { break } }
  }#end x loop
  #Saves final stuff
  temp.align<-lapply(x3, FUN = function(x) paste(x, collapse = ""))
  align.out<-DNAStringSet(unlist(temp.align))
  names(align.out)<-save.names
  return(align.out)
}

slice.trim<-function(input.align, slice.size.bp = 100, threshold = 0.45){

  #makes consensus sequence for comparison
  #input.align<-trimal.align
  input.con<-make.consensus(input.align, method = "majority")
  names(input.con)<-"Reference_Locus"

  comb.align<-append(input.align, input.con)

  #Gets slice information ready
  slice.no<-ceiling(max(width(input.align))/slice.size.bp)
  slice.start<-1
  slice.end<-slice.size.bp

  #checks to see if its out of bounds
  if (slice.end > max(width(input.align))){
    slice.end<-max(width(input.align))
  }#end if check
  output.align<-DNAStringSet()
  for (x in 1:slice.no){

    #Slice alignment into number of slices
    sliced.align<-subseq(comb.align, start = slice.start, end = slice.end)
    #Checks for badly aligned sequences
    bad.align<-pairwise.inf.sites(sliced.align, "Reference_Locus")
    #Remove bad sequence chunks
    rem.seqs<-bad.align[bad.align >= threshold]
    good.align<-sliced.align[!names(sliced.align) %in% names(rem.seqs)]
    #Makes replacement gap seqs for the bad ones
    blank.align<-DNAStringSet()
    if (length(rem.seqs) != 0){
      for (y in 1:length(rem.seqs)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", slice.end-slice.start+1), collapse = "")) )
      }
      names(blank.align)<-names(rem.seqs)
    }#end rem seqs if

    #Saves the slices and cats
    save.slice<-append(good.align, blank.align)
    save.slice<-save.slice[order(names(save.slice))]
    save.names<-names(save.slice)
    output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
    names(output.align)<-save.names

    #Gets new start and stop
    slice.start<-slice.start+100
    slice.end<-slice.end+100
    #checks to see if its out of bounds
    if (slice.end > max(width(input.align))){
      slice.end<-max(width(input.align))
      if (slice.end-slice.start <= 25){ break } else {
        save.slice<-subseq(comb.align, start = slice.start, end = slice.end)
        save.slice<-save.slice[order(names(save.slice))]
        save.names<-names(save.slice)
        output.align<-DNAStringSet(paste0(as.character(output.align), as.character(save.slice)))
        names(output.align)<-save.names
        break
      }
    }#end if
  }#end x loop

  #Removes reference
  output.align<-output.align[names(output.align) != "Reference_Locus"]
  #removes gap only taxa
  str.splitted<-strsplit(as.character(output.align), "")
  x.align<-as.matrix(as.DNAbin(str.splitted) )
  len.temp<-as.character(as.list(x.align))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= 20]
  return.align<-output.align[!names(output.align) %in% names(spp.rem)]
  return(return.align)
}#end FUNCTION

run.trimal<-function(input.align, method = "auto"){
  # input.align<-rem.align
  #Finds probes that match to two or more contigs
  save.rownames<-names(input.align)
  write.align<-as.list(as.character(input.align))
  write.fasta(sequences = write.align, names = names(write.align),
              paste0(gsub(pattern = "\\..*", "", locus.names[i]), ".fa"), nbchar = 1000000, as.string = T)

  input.file<-paste0(gsub(pattern = "\\..*", "", locus.names[i]), ".fa")

  system(paste0(trimal.path, " -in ", input.file, " -out ", input.file, "-tm -automated1"))
  system(paste0("rm ", input.file))

  if (file.exists(paste0(input.file, "-tm")) == F) {
    system(paste0("rm ", input.file))
    print(paste0("deleted. Not enough overlapping data in alignment.") )
    return(DNAStringSet())
  } else { system(paste0("mv ", input.file, "-tm ", input.file)) }

  out.align<-scanFa(FaFile(input.file))

  #Fixes any terrible NA names introduced by trimal
  new.names<-c()
  for (j in 1:length(names(out.align))){
    new.names[j]<-save.rownames[grep(pattern = paste0(names(out.align)[j], "$"), x = save.rownames)]
  }

  temp<-names(out.align)[is.na(names(out.align)) == T]
  if (length(temp) > 0){ stop("there are NAs in the names") }
  names(out.align)<-new.names
  unlink(input.file)
  return(out.align)
}#end function

make.consensus<- function (input.alignment, method = c("majority", "threshold", "IUPAC",
                                                       "profile"), threshold = 0.6, warn.non.IUPAC = FALSE, type = c("DNA", "RNA")) {

  #input.alignment<-trimmed
  #Converts alignment to matrix of characters to be used
  new.align<-strsplit(as.character(input.alignment), "")
  align.in<-matrix(unlist(new.align), ncol = length(new.align[[1]]), byrow = T)

  #Does based on method
  method <- match.arg(method)

  if (method == "IUPAC") {
    type <- match.arg(type)
    res <- apply(align.in, 2, bma, warn.non.IUPAC = warn.non.IUPAC,
                 type = type)
    names(res) <- NULL
  }
  if (method == "majority") {
    majority <- function(x) names(which.max(table(x)))
    res <- apply(align.in, 2, majority)
    names(res) <- NULL
  }
  if (method == "profile") {
    obsvalue <- levels(factor(align.in))
    nrow <- length(obsvalue)
    row.names(align.in) <- NULL
    res <- apply(matali, 2, function(x) table(factor(x, levels = obsvalue)))
  }
  if (method == "threshold") {
    profile <- consensus(align.in, method = "profile")
    profile.rf <- apply(profile, 2, function(x) x/sum(x))
    res <- rownames(profile.rf)[apply(profile.rf, 2, which.max)]
    res <- ifelse(apply(profile.rf, 2, max) >= threshold,
                  res, NA)
    names(res) <- NULL
  }

  out.consensus<-DNAStringSet(paste0(res, collapse = ""))
  names(out.consensus)<-"Consensus_Sequence"
  return(out.consensus)
}

################################################################################################
################################################################################################
############## Step 1. Group exons + introns all markers together   ############################
################################################################################################
################################################################################################
################################################################################################

#Sets up new directory for this stuff
trim.dir = "gene-complete_trimmed"
all.dir = "all-markers_trimmed"
dir.create(paste0(align.dir, "/", trim.dir))
setwd(paste0(align.dir, "/", all.dir))

#Finds loci from the file names
all.names = list.files(paste0(align.dir, "/", all.dir))
file.names = list.files(paste0(align.dir, "/", all.dir))
temp.names = gsub("_exon.*", "", file.names)
temp.names = gsub("-ex.*", "", temp.names)
locus.names = unique(temp.names[duplicated(temp.names)])
gene.names = locus.names

# #Fixes broken loci
# if (length(locus.names) <= 10){
#   exon.names<-file.names[grep("-ex", file.names)]
#   temp.names<-str_split(exon.names, "-")
#   locus.names <- unique(unlist(lapply(temp.names, "[[", 3)))
# }

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){

  ##############
  #STEP 1: Get the taxa in the alignments for these loci
  ##############
  #Current probe sets with hyphens
  if (length(grep(".phy$", locus.names[i])) >= 1){ next }
  setwd(paste0(align.dir, "/", all.dir))
  locus.files = all.names[grep(paste0(locus.names[i]), all.names)]
  locus.files = locus.files[order(gsub(".*-ex", "", locus.files))]

  #Older probe sets with underscore
  if (length(locus.files) == 0){
    locus.files = all.names[grep(paste0("_", locus.names[i], "-"), all.names)]
    locus.files = locus.files[order(gsub(".*-ex", "", locus.files))]
  } #end if

  taxa.names<-c()
  for (j in 1:length(locus.files)){
    #Reads in files
    align<-readAAMultipleAlignment(file = paste(align.dir, "/", all.dir, "/", locus.files[j], sep =""), format = "phylip")
    taxa.names<-append(taxa.names, rownames(align))
  }

  taxa.names<-unique(taxa.names)

  ##############
  #STEP 2: Get different loci
  ##############

  combined.align<-DNAStringSet()
  for (j in 1:length(locus.files)){
    #Reads in files
    align<-readAAMultipleAlignment(file = paste0(align.dir, "/", all.dir, "/", locus.files[j]), format = "phylip")
    align<-DNAStringSet(align)

    add.taxa<-taxa.names[!taxa.names %in% names(align)]

    blank.align<-DNAStringSet()
    if (length(add.taxa) != 0){
      for (y in 1:length(add.taxa)){
        blank.align<-append(blank.align, DNAStringSet(paste0(rep("-", max(width(align))), collapse = "")) )
      }
      names(blank.align)<-add.taxa
    }#end rem seqs if

    #Saves the slices and cats
    new.align<-append(align, blank.align)
    new.align<-new.align[order(names(new.align))]
    save.names<-names(new.align)
    combined.align<-DNAStringSet(paste0(as.character(combined.align), as.character(new.align)))
    names(combined.align)<-save.names
  }#end j loop

  ##############
  #STEP 3: Cleanup and save
  ##############

  #string splitting
  write.temp<-strsplit(as.character(combined.align), "")
  aligned.set<-as.matrix(as.DNAbin(write.temp) )

  #Removes samples that too short individually
  len.temp<-as.character(as.list(aligned.set))
  len.loci<-lapply(len.temp, function (x) x[x != "-"])
  spp.len<-unlist(lapply(len.loci, function (x) length(x)))
  spp.rem<-spp.len[spp.len <= max(spp.len) * as.numeric(min.cov)]
  if (length(spp.rem) > 0){  aligned.set<-aligned.set[!rownames(aligned.set) %in% unique(names(spp.rem)),] }

  #removes loci with too few taxa
  if (length(rownames(aligned.set)) <= as.numeric(min.taxa)){
    print(paste(locus.names[i], " deleted. Too few taxa after trimming.", sep = ""))
    next
  }

  #readies for saving
  write.phy(aligned.set, file= paste0(align.dir, "/", trim.dir, "/", locus.names[i], ".phy"), interleave = F)

}#end i loop


#Sets up new directory for this stuff
new.trim.dir = "gene-all-markers_trimmed"
setwd(paste0(align.dir))
all.dir = "all-markers_trimmed"
system(paste0("cp -r ", align.dir, "/", trim.dir, " ",
              align.dir, "/", new.trim.dir))

#Finds loci from the file names
all.names = list.files(paste0(align.dir, "/", all.dir))
file.names = list.files(paste0(align.dir, "/", all.dir))
temp.names = gsub("_exon.*", "", file.names)
temp.names = gsub("-ex.*", "", temp.names)
locus.names = temp.names[!temp.names %in% gene.names]

#Loops through each locus and does operations on them
for (i in 1:length(locus.names)){

  copy.file = file.names[grep(locus.names[i], file.names)]

  if (length(copy.file) != 1){stop("too mayn files")}

  system(paste0("cp ", align.dir, "/", all.dir, "/", copy.file, " ",
                align.dir, "/", new.trim.dir, "/", copy.file))


}#end i loop






# END SCRIPT

