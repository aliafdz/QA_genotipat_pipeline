
library(stringr)
library(Biostrings)
source("./R/seqanalfns.v4.5.R")

###  Fitxers a tractar i descripció
#####################################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
o <- order(FlTbl$Pat.ID,FlTbl$Ampl.Nm,FlTbl$Str)
FlTbl <- FlTbl[o,]
###  Compose fasta file names
idx.fw <- which(FlTbl$Str=="fw")
idx.rv <- which(FlTbl$Str=="rv")

flnms <- paste("MACHpl02",FlTbl$Pat.ID[idx.fw],
               FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
flnms <- file.path(mach.Dir,flnms)
pnms <- paste(FlTbl$Pat.ID[idx.fw],FlTbl$Ampl.Nm[idx.fw],sep=" - ")

pdf(file=file.path(repDir,"GapsBarPlots.pdf"),paper="a4",width=7,
    height=10)
par(mfrow=c(4,1))
par(mar=c(4.5, 4, 3, 4) + 0.1)

for(i in 1:length(flnms))
{ 
  if(!file.exists(flnms[i])) next
  lst <- read.ampl.seqs(flnms[i],mnr=1)
  imast <- which.max(lst$IDs$nseqs)
  rsq <- lst$seqs[imast]
  xmx <- nchar(rsq)
  seqs <- lst$seqs
  nr <- lst$IDs$nseqs
  
  rnt <- strsplit(rsq,split="")[[1]]  
  tp <- integer(length(rnt))
  ins.pos <- rnt=="-"
  tp[ins.pos] <- 1  # insercions

  ntmat <- t(sapply(seqs,function(s) strsplit(s,split="")[[1]]))
  del.pos <- apply(ntmat,2,function(vnt) any(vnt=="-"))
  tp[!ins.pos & del.pos] <- 2  # deleccions
  tp <- tp+1
  
  nb <- integer(length(rnt))
  if( sum(ins.pos) )
  { nb[ins.pos] <- apply(ntmat[,ins.pos,drop=FALSE],2,function(vnt)
                           sum(nr[vnt!="-"]))
  }
  if( sum(del.pos) )  
  { nb[!ins.pos & del.pos] <- apply(ntmat[,!ins.pos & del.pos,drop=FALSE],2,
                 function(vnt) sum(nr[vnt=="-"]))
  }				 
  plot(nb,type="h",col=c("black","red","blue")[tp],lwd=2,yaxt="n",
       xlab="MA position",ylab="reads",xlim=c(0,xmx),ylim=c(0,max(nb[1:xmx])))
  ytk <- axTicks(2)
  axis(side=2,at=ytk,las=2)
  pct <- round(ytk/sum(nr)*100,2)
  axis(side=4,at=ytk,labels=pct,las=2,cex.axis=0.8)
  mtext("Percentage",side=4,line=3,cex=0.6)
  title(main=pnms[i],line=1)
}

dev.off()
