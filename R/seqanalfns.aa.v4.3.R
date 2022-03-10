
############################################################
###    AMINOACID SEQUENCE ANALYSIS, UTILITY FUNCTIONS    ###
###           source("seqanalfns.aa.v3.5.txt")           ###
############################################################

library(Biostrings)
library(RColorBrewer)
cls <- brewer.pal(8,"Dark2")

###  Aminoacid codes
aa3.nms <- c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His",
         "Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp",
         "Tyr","Val","Xaa","Stp")
aa.nms <- c("A","R","N","D","C","Q","E","G","H","I","L","K",
            "M","F","P","S","T","W","Y","V","X","*")
naas <- length(aa.nms)

###  Aminoacid types and corresponding colors
aa.tp <-  c( 5,  1,  3,  2,  4,  3,  2,  4,  1,  5,  5,  1,
             5,  5,  4,  3,  3,  5,  5,  5,  6,  7)
names(aa.tp) <- aa.nms
aatplv <- c("Positive","Negative","Polar","Special","Hydrophobic",
            "Unknown","STOP")
aatpnm <- c("Pos","Neg","Pol","Spc","Phob","Unk","STP")
aa.cls <- c(cls[1:4],cls[7:8],"Black")

	
###  Table of aa frequencies at position i
aa.PosTbl <- function(i,seqs,w)
{ res <- integer(naas)
  names(res) <- aa.nms
  tbl <- table(substr(seqs,i,i))
  res[names(tbl)] <- tbl
  res
}

###  Table of population aa frequencies at position i
###    given sequence frequencies
aa.PosTbl.w <- function(i,seqs,w)
{ res <- integer(naas)
  names(res) <- aa.nms
  lets <- substr(seqs,i,i)
  for(i in 1:naas)
    res[i] <- sum(ifelse(lets==aa.nms[i],1,0)*w)
  res
}


###  Matrix of aa frequencies per position
aa.SeqsTbl <- function(seqs)
{ tbl <- t(sapply(1:nchar(seqs[1]),function(i) aa.PosTbl(i,seqs)))
  rownames(tbl) <- 1:nchar(seqs[1])
  tbl
}


###  Matrix of population aa frequencies per position
###    given sequence frequencies
aa.SeqsTbl.w <- function(seqs,w)
{ tbl <- t(sapply(1:nchar(seqs[1]),function(i) aa.PosTbl.w(i,seqs,w)))
  rownames(tbl) <- 1:nchar(seqs[1])
  tbl
}


###  Matrix of mutations frequencies per position
###  (the most frequent is set to 0 leaving only the mutations)
aa.MutsTbl <- function(seq.tbl)
{ j <- apply(seq.tbl,1,which.max)
  seq.tbl[cbind(1:nrow(seq.tbl),j)] <- 0
  seq.tbl
}


###  Consensus sequence build from the matrix of aminoacid 
###    frequencies per position. The most frequent aa at 
###    each position is taken.
aa.ConsSeq <- function(seq.tbl)
{ j <- apply(seq.tbl,1,which.max)
  paste(colnames(seq.tbl)[j],collapse="")
}


###  Matrix of summary mutations
aa.SummaryMuts <- function(seqs,off=0)
{ pos.tbl <- aa.SeqsTbl(seqs)
  mut.tbl <- aa.MutsTbl(pos.tbl)
  flags <- apply(mut.tbl,1,sum)>0
  pos <- which(flags)
  res <- cbind(pos=pos+off,pos.tbl[flags,],
               round(pos.tbl[flags,]/sum(pos.tbl[1,])*100,3))
  colnames(res) <- c("pos",paste("f",aa.nms,sep= ""),
                     paste("p",aa.nms,sep= ""))
  rownames(res) <- 1:nrow(res)
  res
}


###  Matrix of summary population mutations
###    given sequence frequencies
aa.SummaryMuts.w <- function(seqs,w,off=0)
{ 
  pos.tbl <- aa.SeqsTbl.w(seqs,w)
  mut.tbl <- aa.MutsTbl(pos.tbl)
  flags <- apply(mut.tbl,1,sum)>0
  pos <- which(flags)
  if(length(pos)>1)
  { res <- cbind(pos=pos+off,pos.tbl[flags,],
               round(pos.tbl[flags,]/sum(pos.tbl[1,])*100,3))
    colnames(res) <- c("pos",paste("f",aa.nms,sep= ""),
                     paste("p",aa.nms,sep= ""))
  } else {
    res <- c(pos+off,pos.tbl[flags,],
               round(pos.tbl[flags,]/sum(pos.tbl[1,])*100,3))
    names(res) <- c("pos",paste("f",aa.nms,sep= ""),
                     paste("p",aa.nms,sep= ""))
    res <- t(data.frame(res))
  }	
  rownames(res) <- 1:nrow(res)
  res
}


###  Segregating sites: Number of sites with mutations
###    (independent of sequence weights)
aa.SegSites <- function(seqs)
{ nrow(aa.SummaryMuts(seqs)) }


###  Segregating sites: Number of sites with mutations
###    given sequence frequencies
aa.SegSites.w <- function(seqs,w)
{ nrow(aa.SummaryMuts.w(seqs)) }


### Total number of mutations (Eta)
###  (independent of sequence weights)
aa.TotalMutations <- function(mut.tbl)
  sum(apply(mut.tbl,1,function(x) sum(x>0)))


### Total number of population mutations
aa.TotalMutations.w <- function(mut.tbl)
  sum( apply(mut.tbl,1,function(x) sum(x)) )


### Number of polymorphic sites (segregating sites) (S)
###  (independent of sequence weights)
aa.PolymorphicSites <- function(mut.tbl)
  sum( apply(mut.tbl,1,sum) > 0 )


###  Aminoacid differences between two sequences (no gaps) 
aa.PairDiffs <- function(seq1,seq2)
  sum( strsplit(seq1,split="")[[1]] != strsplit(seq2,split="")[[1]] )


###  Matrix of pairwise differences
aa.PairwiseDiffs <- function(seqs)
{ n <- length(seqs)
  d <- matrix(0,n,n)
  for(i in 1:(n-1))
    for(j in (i+1):n)
       d[i,j] <- aa.PairDiffs(seqs[i],seqs[j])
  d+t(d)
}


###  Matrix of population pairwise differences
###    given sequence frequencies
aa.PairwiseDiffs.w <- function(seqs,w)
{ n <- length(seqs)
  d <- matrix(0,n,n)
  for(i in 1:(n-1))
    for(j in (i+1):n)
       d[i,j] <- aa.PairDiffs(seqs[i],seqs[j])*w[i]*w[j]
  d+t(d)
}

###  Mean number of aminoacid differences
MeanAaDiffs <- function(d)
{ n <- nrow(d)
  k <- 0
  for(i in 1:(n-1))
    for(j in (i+1):n)
       k <- k + d[i,j]
  k*2/(n*(n-1))
}
 

###  Mean number of population aminoacid differences
###    given sequence frequencies
MeanAaDiffs.w <- function(d,w)
{ n <- nrow(d)
  nseqs <- sum(w)
  k <- 0
  for(i in 1:(n-1))
    for(j in (i+1):n)
       k <- k + d[i,j]
  k*2/(nseqs*(nseqs-1))
}


###  Information content at a site
###    given the aminoacids frequencies
aa.InfContent <- function(v)
{ v <- v/sum(v)
  lgv <- ifelse(v==0,0,log2(v))
  log2(naas-1)+sum(v*lgv)
}


###  Information content at a site given the absolute mutations
###    frequency, and the total nunber of population sequences
aa.InfContent.mut <- function(v,tseq)
{ v <- c(v/tseq,(tseq-sum(v))/tseq)
  lgv <- ifelse(v==0,0,log2(v))
  log2(naas-1)+sum(v*lgv)
}


###  Information content at a site
###    given the mutations percentages
aa.InfContent.pmut <- function(v)
{ v <- v/100
  v <- c(v,1-sum(v))
  lgv <- ifelse(v==0,0,log2(v))
  log2(naas-1)+sum(v*lgv)
}


###  Quasispecies aa Shannon entropy
aa.QSEntropy <- function(w)
{ w <- w/sum(w)
  lgw <- ifelse(w==0,0,log(w))
  -sum(w*lgw)
}


###  Read amplicon aligned sequences 
######################################  
read.ampl.seqs <- function(flnm,mnr=2,mpct=NULL)
{
  seqs <- as.character(readAAStringSet(flnm))
  IDstr <- names(seqs)
  n <- length(IDstr)

  nms <- character(length=n)
  sts <- matrix(0,nrow=n,ncol=2)
  colnames(sts) <- c("nseqs","pct1")
  for(j in 1:n)
  { strs <- strsplit(IDstr[j],split="\\|")[[1]]
    nms[j] <- strs[1]
    sts[j,] <- as.numeric(strs[2:3])
  }
  IDs <- data.frame(ID=nms,sts,stringsAsFactors=FALSE)
  nall <- nrow(IDs)
  tnr <- sum(IDs$nseqs)
 
  ###  Filter by minimum reads by haplotype
  if(is.null(mpct))
  { flags <- IDs$nseqs >= mnr
  } else {
    flags <- IDs$pct1 >= mpct & IDs$nseqs >= mnr
  }
  return(list(IDs=IDs[flags,],seqs=seqs[flags],nall=nall,tnr=tnr))
}


###  Print (in short) the aa mutations observed in each position
print.aa.polym <- function(mtbl)
{ 
  if(is.null(dim(mtbl)))
  { x <- mtbl[(1:naas)+1]
    names(x) <- 1:naas
    x <- sort(x,decreasing=TRUE)
    ks <- as.integer(names(x))
    cat("\n",mtbl[1],": ",sep="")
    cat( sapply(1:naas,function(a)
            ifelse(x[a]>0,
              paste(aa3.nms[ks[a]]," ",x[a],", ",sep=""), "")),
         sep="")
    cat("\n")
    return()
  }
  for(i in 1:nrow(mtbl))
  { x <- mtbl[i,(1:naas)+1]
    names(x) <- 1:naas
    x <- sort(x,decreasing=TRUE)
    ks <- as.integer(names(x))
    cat("\n",mtbl[i,1],": ",sep="")
    cat( sapply(1:naas,function(a)
            ifelse(x[a]>0,
              paste(aa3.nms[ks[a]]," ",x[a],", ",sep=""), "")),
         sep="")
  }
  cat("\n")
}


###  Read aa amplicon aligned sequences 
##########################################  
read.aa.ampl.seqs <- function(flnm,mnr=2)
{
  seqs <- as.character(readAAStringSet(flnm))
  IDstr <- names(seqs)
  n <- length(IDstr)

  nms <- character(length=n)
  sts <- matrix(0,nrow=n,ncol=2)
  colnames(sts) <- c("nseqs","pct1")
  for(j in 1:n)
  { strs <- strsplit(IDstr[j],split="\\|")[[1]]
    nms[j] <- strs[1]
    sts[j,] <- as.numeric(strs[2:3])
  }
  IDs <- data.frame(ID=nms,sts,stringsAsFactors=FALSE)
  nall <- nrow(IDs)
  tnr <- sum(IDs$nseqs)
 
  ###  Filter by minimum reads by haplotype
  flags <- IDs$nseqs >= mnr
  return(list(IDs=IDs[flags,],seqs=seqs[flags],nall=nall,tnr=tnr))
}


### After right and/or left trimming recollapse identical sequences
##################################################################### 
CollapseSeqs <- function(seqs,IDs)
{
  for(i in 1:(length(seqs)-1))
  { if(IDs$nseqs[i]==0) next
    for(j in (i+1):length(seqs))
    { if(IDs$nseqs[i] == 0 | IDs$nseqs[j] == 0) next
      if(seqs[i]==seqs[j]) 
      { IDs$nseqs[i] <- IDs$nseqs[i]+IDs$nseqs[j]
        IDs$nseqs[j] <- 0 
      } 
    }
  }
  flags <- IDs$nseqs > 0
  list(seqs=seqs[flags],IDs=IDs[flags,])
}


###  Print formated sequence
##############################
print.seq <- function(rseq)
{
  ln <- nchar(rseq)
  lt <- strsplit(rseq,split="")[[1]]
  nr <- ln%/%50
  if(ln%%50) nr <- nr+1
  u <- paste(rep("1234567890",5),collapse="")
  d1 <- "    .    1    .    2    .    3    .    4    .    5"
  d2 <- "    .    6    .    7    .    8    .    9    .    0"
  k <- 1
  j <- 1
  for(i in 1:nr)
  { lst <- k+50-1
    if(lst>ln) lst <- ln
    str <- paste(lt[k:lst],collapse="")
    strln <- lst-k+1
    d <- ifelse(i%%2,d1,d2)
    cat("\n",substr(d,1,strln))
    cat("\n",substr(u,1,strln))
    cat("\n",str)
    k <- lst+1
  }
  cat("\n\n")
}

           ###   END OF seqanalfns.aa.v4.3.R   ###
