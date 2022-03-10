############################################################
###    RAW CONSENSUS HAPLOTYPES by MULTIPLE ALIGNMENT    ###
############################################################


###  Formatar sencers emplenant amb 0 per l'esquerra
######################################################
zeroFillInt2Char <- function(x,ln)
{ x <- paste("000000",x,sep="")
  substring(x,nchar(x)-ln+1)
}


###  Split seq names with counts and percentage 
#################################################
split.fasta.names <- function(seqs,mnr=2)
{ IDstr <- names(seqs)
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


###  Read amplicon aligned sequences 
#######################################
read.ampl.seqs <- function(flnm,mnr=2)
{ seqs <- as.character(readDNAStringSet(flnm))
  return( split.fasta.names(seqs,mnr) )
}


###  Escriu unes seqüencies a un fitxer fasta
###############################################
write.fasta <- function(seqs,flnm)
{ writeXStringSet(seqs,flnm) }


###  Llegeix les seqüencies d'un fitxer fasta
###############################################
read.fasta <- function(flnm)
{ readDNAStringSet(flnm) }


###  Aliniament múltiple per muscle de les seqs
###  Les seqüècies s'entren i es tornen com a DNAStringSet
############################################################
#muscle <- "C:\\Muscle\\muscle3.8.31_i86win32.exe"
#muscle <- "D:\\UltraSeq\\muscle\\muscle3.8.31_i86win32.exe"
#muscle <- "E:\\UltraSeq\\muscle\\muscle3.8.31_i86win32.exe"
#muscle <- "C:\\Docs\\UltraSeq\\muscle\\muscle3.8.31_i86win32.exe"
muscle.cl.opts <- c("-log muscle.log")
doMuscle <- function(seqs)
{ tmp.file <- file.path(tmp.Dir,"muscleInFile.fna")
  res.file <- file.path(tmp.Dir,"muscleOutFile.fna")
  if(file.exists(res.file)) file.remove(res.file)
  write.fasta( seqs, tmp.file )
  in.file <- paste("-in ",tmp.file,sep="")
  out.file <- paste("-out ",res.file,sep="")
  command <- paste(muscle,in.file,out.file,
                paste(muscle.cl.opts,collapse=" "),sep=" ")
  system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
          ignore.stdout=FALSE,invisible=TRUE)
  if( file.exists(res.file) )
    return( read.fasta(res.file) )
  return(NULL)
}


###  Taula de freqüències de nucleòtid per posició
#####################################################
SeqsPosTable <- function(seqs,nr)
{ col.freqs<- function(s,nr) 
  { c( sum(nr[s=="A"]),sum(nr[s=="C"]),sum(nr[s=="G"]),
       sum(nr[s=="T"]),sum(nr[s=="-"]) )
  }
  nuc.mat <- t(sapply(seqs,function(s) strsplit(s,split="")[[1]]))
  frpos <- t(apply(nuc.mat,2,col.freqs,nr))
  colnames(frpos) <- c("A","C","G","T","-")
  frpos
}

       
###  Alinear les distribucions d'haplotips de dues poblacions
###############################################################
PopsAlgnHist <- function(nA,seqsA,nB,seqsB)
{ names(nA) <- seqsA
  names(nB) <- seqsB
  nms <- union(seqsA,seqsB)    # tots: A + B
  nb <- length(nms)
  pA <- integer(nb)
  idx <- which(nms %in% seqsA) # els de A en tots
  pA[idx] <- nA[nms[idx]]  
  pB <- integer(nb)
  idx <- which(nms %in% seqsB) # els de B en tots
  pB[idx] <- nB[nms[idx]]  
  list(pA=pA,pB=pB,Hpl=nms)
}


###  Confrontar haplotips a FW i RV
Intersect.FWRV <- function(nA,seqsA,nB,seqsB)  
{
  ###  Aliniar distribucions
  lst <- PopsAlgnHist(nA,seqsA,nB,seqsB)
  ###  Haplotips comuns i solapament global
  fl <- lst$pA>0 & lst$pB>0
  ov.a <- sum(lst$pA[fl]+lst$pB[fl])/sum(lst$pA+lst$pB)
  ###  Frequencies i solapament per intersecció
  pFW <- (lst$pA/sum(lst$pA))[fl]
  pRV <- (lst$pB/sum(lst$pB))[fl]
  p <- pmin(pFW,pRV) # interseccio
  ov.i <- sum(p)  
  p <- p/sum(p)      # renormalitzacio

  list(p=p,seqs=lst$Hpl[fl],ov.i=ov.i,ov.a=ov.a,pA=lst$pA,pB=lst$pB)
}


### Representa les distribucions aliniades dels haplotips
###########################################################
PlotHplHistos <- function(tt,pA,pB,p)
{ pA <- pA/sum(pA)
  pB <- pB/sum(pB)
  ymx <- max(c(pA,pB))
  barplot(pA,ylim=c(0,ymx)); abline(h=0)
  title(main="FW strand haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  barplot(pB,ylim=c(0,ymx)); abline(h=0)
  title(main="RV strand haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  barplot(p); abline(h=0)
  title(main="Intersected haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
}



###  Dóna nom als haplotipus en funció del nombre de mutacions respecte
###   la màster i de la seva frequencia poblacional, i salva les 
###   sequencies en un fitxer fasta.
########################################################################
SaveHaplotypes <- function(flnm,bseqs,nr,max.difs=250)  #,seq0)
{
  ##  Si només hi ha una seq.
  if(length(bseqs)==1)
  { bnts <- strsplit(bseqs[1],split="")[[1]]
	bnts <- bnts[bnts!="-"]
	bseqs <- paste(bnts,collapse="")
	names(bseqs) <- paste("Hpl.0.0001",nr,100,sep="|")
    write.fasta(DNAStringSet(bseqs),flnm)
	return( list(bseqs=bseqs,nr=nr,nm=0) )
  }

  ##  Determinar diferències respecte la màster
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[which.max(nr)])
  nm <- nmismatch(psa)

  ##  Eliminar seqs amb massa diferències
  bseqs <- bseqs[nm<=max.difs]
  nr <- nr[nm<=max.difs]
  nm <- nm[nm<=max.difs]
  
  if(length(bseqs)==1)
  { bnts <- strsplit(bseqs[1],split="")[[1]]
	bnts <- bnts[bnts!="-"]
	bseqs <- paste(bnts,collapse="")
	names(bseqs) <- paste("Hpl.0.0001",nr,100,sep="|")
    write.fasta(DNAStringSet(bseqs),flnm)
	return( list(bseqs=bseqs,nr=nr,nm=0) )
  }

  ##  Ordenar per nombre de diferències
  o <- order(nm)
  bseqs <- bseqs[o]
  nr <- nr[o]
  nm <- nm[o]

  ##  Numero d'ordre dins de cada nombre de mutacions
  tnm <- table(nm)
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  
  ##  Ordenar per frequencia descendent dins de cada nombre de mutacions
  for(i in as.integer(names(tnm)))
  { idx <- which(nm==i)
    o <- order(nr[idx],decreasing=TRUE)
    bseqs[idx] <- bseqs[idx[o]]
    nr[idx] <- nr[idx[o]]
  }

  ##  Calcular frequencia relativa
  frq <- round(nr/sum(nr)*100,2)

  ## Nom complet per cada haplotipus
  nms <- paste("Hpl",nm,zeroFillInt2Char(isq,4),sep=".")

  ##  Capçalera fasta amb nom, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ###  Afegir-hr la RefSeq
  #bseqs <- c(seq0,bseqs)

  ###  Eliminar columnes de tot gaps
  nuc.mat <- t(sapply(bseqs,function(s) strsplit(s,split="")[[1]]))
  fl <- apply(nuc.mat,2,function(st) all(st=="-"))
  if(sum(fl))
  { nuc.mat <- nuc.mat[,!fl]
    bseqs <- apply(nuc.mat,1,paste,collapse="")
  }

  ##  Salvar a fasta
  write.fasta(DNAStringSet(bseqs),flnm)

  list(bseqs=bseqs,nr=nr,nm=nm)
}


###  Mesues de diversitat d'un amplicó
########################################
AmplStats <- function(seqs,nr)
{
  smres <- data.frame( Hapl=integer(1),Eta=integer(1),
                       S=integer(1),Mf=numeric(1),Sn=numeric(1),
                       Pi=numeric(1) )
  smres$Hapl <- length(seqs)
  ###  Compute the nucleotide frequencies and the mutations table
  pos.tbl <- SeqsTbl.w(seqs,nr)
  mut.tbl <- MutsTbl(pos.tbl)
  ### Total number of mutations (Eta)
  smres$Eta <- TotalMutations(mut.tbl)
  ### Number of polymorphic sites (segregating sites, S)
  smres$S <- PolymorphicSites(mut.tbl)
  ### Mutations frequency
  smres$Mf <- TotalMutations.w(mut.tbl) / 
                (sum(nr)*nchar(seqs[1]))
  ### Normalized quasispecies Shannon entropy
  smres$Sn <- QSEntropy(nr)/log(smres$Hapl)
  ###  Pairwise differences
  d.ij <- PairwiseDiffs.w(seqs,nr)
  ###  Nucleotide diversity (Pi)
  smres$Pi <- MeanNuclDiffs.w(d.ij,nr)/nchar(seqs[1])
  smres
}

#-------------------------------------------------------------------#

library(Biostrings)
library(ape)
source("./R/seqanalfns.v4.5.R")

###  Llegir la descripció de mostres
#######################################
samples <- read.table(file.path(desc.Dir,"samples.csv"), sep="\t",
                      header=T,stringsAsFactors=F)
primers <- read.table(file.path(desc.Dir,"primers.csv"), sep="\t",
                      header=T,stringsAsFactors=F)
# RefSeqs <- read.table(file.path(desc.Dir,"RefSeqs.csv"), sep="\t",
                      # header=T,stringsAsFactor=F)

#  Fitxers a tractar
########################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
o <- order(FlTbl$Pat.ID,FlTbl$Ampl.Nm,FlTbl$Str)
FlTbl <- FlTbl[o,]
###  Compose fasta file names
in.files <- file.path(trimDir,FlTbl$File.Name)
idx.fw <- which(FlTbl$Str=="fw")
idx.rv <- which(FlTbl$Str=="rv")

out.flnms <- paste("MACHpl02",FlTbl$Pat.ID[idx.fw],
               FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
out.flnms <- file.path(mach.Dir,out.flnms)

ma.flnms <- paste("MAfwrv",FlTbl$Pat.ID[idx.fw],
               FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
ma.flnms <- file.path(mach.Dir,ma.flnms)

n <- length(in.files)/2

rdf.fw <- data.frame(fw.all=integer(n),fw.lowf=integer(n),
                     fw.in=integer(n),fw.unq=integer(n),fw.com=integer(n))
rdf.rv <- data.frame(rv.all=integer(n),rv.lowf=integer(n),
                     rv.in=integer(n),rv.unq=integer(n),rv.com=integer(n))
rdf.gbl <- data.frame(all=integer(n),lowf=integer(n),unq=integer(n),
                      ovrlp=numeric(n),common=numeric(n),Fn.rd=integer(n))

pdf(file.path(repDir,"MA.Intersects.plots.pdf"),paper="a4",
    width=7,height=11)
par(mfrow=c(3,1))

for(i in 1:n)
{ 
  if(!file.exists(in.files[idx.fw[i]])) next
  if(!file.exists(in.files[idx.rv[i]])) next

  # Filtrat implicit per mnr = min.rd
  lst1 <- read.ampl.seqs(in.files[idx.fw[i]],mnr=min.rd)
  nr1 <- lst1$IDs$nseqs
  seqs <- lst1$seqs
  # Eliminar les curtes
  fl <- nchar(seqs) >= min.seq.len
  nr1 <- nr1[fl]
  seqs <- seqs[fl]
  # Filtrar per mínima abundància
  fl1 <- nr1/sum(nr1)*100 >= a.cut
  rdf.fw$fw.all[i]  <- sum(nr1)
  rdf.fw$fw.lowf[i] <- sum(nr1[!fl1])
  seqs <- seqs[fl1]
  names(seqs) <- sub("Hpl","HplFw",names(seqs))
  aseqs <- seqs
  rawln <- nchar(seqs[1])

  # Filtrat implicit per mnr = min.rd
  lst2 <- read.ampl.seqs(in.files[idx.rv[i]],mnr=min.rd) 
  nr2 <- lst2$IDs$nseqs
  seqs <- lst2$seqs
  # Eliminar les curtes
  fl <- nchar(seqs) >= min.seq.len
  nr2 <- nr2[fl]
  seqs <- seqs[fl]
  # Filtrar per mínima abundància
  fl2 <- nr2/sum(nr2)*100 >= a.cut
  rdf.rv$rv.all[i]  <- sum(nr2)
  rdf.rv$rv.lowf[i] <- sum(nr2[!fl2])
  seqs <- seqs[fl2]
  names(seqs) <- sub("Hpl","HplRv",names(seqs))
  aseqs <- c(aseqs,seqs)

  ipr <- FlTbl$Pr.ID[idx.fw[i]]

  ###  Aliniament múltiple per muscle de hpl fw, hpl rv i RefSeq
  seqs <- doMuscle(DNAStringSet(aseqs))
  file.copy(file.path(tmp.Dir,"muscleOutFile.fna"),ma.flnms[i],
            overwrite=TRUE)

  lst <- split.fasta.names(as.character(seqs))
  nr <- lst$IDs$nseqs
  seqs <- lst$seqs

  ###  Separar FW i RV ja aliniades
  ifw <- grep("^HplFw",names(seqs))
  irv <- grep("^HplRv",names(seqs))

  lst <- Intersect.FWRV(nr[ifw],seqs[ifw],nr[irv],seqs[irv])
  fl <- lst$pA>0 & lst$pB>0  ###  Haplos comuns

  # cat("\n  Haplotypes distribution overlap: ",round(lst$ov.i*100,2),"%",
      # sep="")
  # cat("\n       Overlap as reads in common: ",round(lst$ov.a*100,2),"%\n",
      # sep="")
  # cat("\n       Input reads: ",sum(lst$pA),", ",sum(lst$pB),sep="")
  # cat("\n   Reads in common: ",sum(lst$pA[fl]),", ",sum(lst$pB[fl]),sep="")
  # cat("\n        Reads lost: ",sum(lst$pA[!fl]),", ",
                               # sum(lst$pB[!fl]),sep="")
  # cat("\n\n")

  rdf.fw$fw.in[i]  <- sum(lst$pA)
  rdf.rv$rv.in[i]  <- sum(lst$pB)
  rdf.fw$fw.com[i] <- sum(lst$pA[fl])
  rdf.rv$rv.com[i] <- sum(lst$pB[fl])
  rdf.fw$fw.unq[i] <- sum(lst$pA[!fl])
  rdf.rv$rv.unq[i] <- sum(lst$pB[!fl])

  rdf.gbl$all[i]    <- rdf.fw$fw.all[i]+rdf.rv$rv.all[i]
  rdf.gbl$ovrlp[i]  <- round(lst$ov.i*100,2)
  rdf.gbl$common[i] <- round(lst$ov.a*100,2)
  rdf.gbl$Fn.rd[i]  <- rdf.fw$fw.com[i]+rdf.rv$rv.com[i]
  rdf.gbl$lowf[i]   <- rdf.fw$fw.lowf[i] + rdf.rv$rv.lowf[i]
  rdf.gbl$unq[i]    <- rdf.fw$fw.unq[i] + rdf.rv$rv.unq[i]

  ###  No overlap, skip.
  if(sum(fl)==0) 
  { cat("\n--------------------------------------------------\n")
    next
  }

  ###  Plot aligned haplotypes bar plots
  tt <- paste(FlTbl$Pat.ID[idx.fw[i]]," - ",FlTbl$Ampl.Nm[idx.fw[i]],sep="")
  p <- lst$pA+lst$pB
  p[ lst$pA==0 | lst$pB==0 ] <- 0
  p <- p/sum(p) 
  PlotHplHistos(tt,lst$pA,lst$pB,p)

  ###  Simple suma de reads a FW i RV, proporcio com a promig ponderat
  rds <- lst$pA[fl]+lst$pB[fl]

  ###  Save common haplos with estimated frequencies
  SaveHaplotypes(out.flnms[i],lst$seqs,rds)  #,seq0)
}

dev.off()


sink(file=file.path(repDir,"MA.Intersects-SummRprt.txt"))

cat("\n   FW + RV HAPLOTYPES INTERSECTIONS")
cat("\n======================================\n")
cat("\nCutting FW and RV at ",a.cut,"% followed by haplotypes intersection.\n",
    sep="")
cat("\nFrequencies as sum of reads FW+RV.\n\n")

cat("\nSUMMARY RESULTS BY READS NUMBER")
cat("\n===============================\n\n")
fl.fw <- FlTbl$Str=="fw"
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.fw,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.rv,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.gbl,stringsAsFactors=FALSE)
print(frdf)

cat("\n\nTotal counts:\n\n")
tots.fw <- apply(rdf.fw,2,sum)
print(tots.fw)
cat("\n")
tots.rv <- apply(rdf.rv,2,sum)
print(tots.rv)
cat("\n")
tots.gbl <- apply(rdf.gbl[,c(1:3,6)],2,sum)
print(tots.gbl)
cat("\n")

cat("\nPercentage:\n\n")
ptts <- round(tots.fw/tots.fw[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.rv/tots.rv[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.gbl/tots.gbl[1]*100,2)
print(ptts)
cat("\n")

sink()

save(frdf,file=file.path(repDir,"IntersectedReads.RData"))

pdf.flnm2 <- "IntersectBarplots.pdf"
pdf(file.path(repDir,pdf.flnm2),paper="a4r",width=10,height=6)
par(mar=c(7,4,4,2)+.1)

pal <- cls[1:2]
bnms <- paste(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw])
vals <- cbind(rdf.fw$fw.com,rdf.rv$rv.com)
ymx <- max(rowSums(vals))*1.2
barplot(t(vals),col=pal,ylim=c(0,ymx),names.arg=bnms,
        las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("fw","rv"),cex=0.8)
title(main="Intersected reads")

vals <- cbind(rowSums(rdf.fw[,c("fw.lowf","fw.unq")]),
              rowSums(rdf.rv[,c("rv.lowf","rv.unq")]))
ymx <- max(rowSums(vals))*1.2
barplot(t(vals),col=pal,ylim=c(0,ymx),names.arg=bnms,
        las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("fw","rv"),cex=0.8)
title(main="Filtered out reads")
			  
pal <- cls[c(1,3,2,4)]
mbp <- data.matrix(frdf[,c(8,5,4)])
mbp.nms <- paste(frdf$Pat.ID,frdf$Ampl.Nm,sep=".")
omar <- par(mar=c(6,3.5,3,2))
barplot(t(mbp),col=pal,ylim=c(0,max(rowSums(mbp))*1.2),names.arg=mbp.nms,
    las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("Common","Unique","Low freq."),
        cex=0.8)
par(mar=omar)
barplot(frdf$common,col="lavender",ylim=c(0,100),names.arg=mbp.nms,
        las=2,cex.names=0.6,cex.axis=0.8)
title(main="FW + RV intersection yield",line=1)

###  Carregar PoolTbl i Fltbl 
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))

all.nms <- paste(FlTbl$Pat.ID,FlTbl$Ampl.Nm,sep=".")
trds <- tapply(FlTbl$Reads,all.nms,sum)
all.nms <- names(trds)
frds <- frdf$Fn.rd
names(frds) <- paste(frdf$Pat.ID,frdf$Ampl.Nm,sep=".")
ytbl <- matrix(0,nrow=length(all.nms),ncol=2)
rownames(ytbl) <- all.nms
colnames(ytbl) <- c("All","Passed")
ytbl[,1] <- trds
ytbl[names(frds),2] <- frds
omar <- par(mar=c(6,3.5,3,2))
barplot(t(ytbl),beside=TRUE,col=cls[1:2],ylim=c(0,max(ytbl)*1.2),las=2,
     names.arg=rownames(ytbl),cex.names=0.6,cex.axis=0.8)
abline(h=0)
legend("top",horiz=TRUE,fill=cls,legend=colnames(ytbl),cex=0.8)
par(mar=omar)
barplot(as.vector(ytbl[,2]/ytbl[,1]*100),col="lavender",ylim=c(0,100),
        names.arg=rownames(ytbl),las=2,cex.names=0.6,cex.axis=0.8)
abline(h=0)
title(main="Global yield",line=1)

dev.off()


# ###  Global yield by step

# load(file=file.path(repDir,"FLASH_table.RData"))
# rownames(flash.res) <- str_extract(rownames(flash.res),"^[A-Za-z0-9\\.-]+")

# p.nms <- rownames(PoolTbl)
# flash.res <- flash.res[p.nms,]
# rds.raw <- flash.res[,1]+flash.res[,2]
# rds.flash <- flash.res[,1]
# rds.MID <- PoolTbl[p.nms,"MIDReads"]
# rds.dmult <- PoolTbl[p.nms,"PrimerReads"]

# frdf.pool <- sapply(1:nrow(frdf), function(i)
               # samples$Pool.Nm[ which(samples$Patient.ID==frdf$Pat.ID[i] &
		                              # samples$Primer.ID==frdf$Ampl.Nm[i])[1] ])
								  
# rds.filt <- tapply(frdf$all,frdf.pool,sum)[p.nms]
# rds.ints <- tapply(frdf$Fn.rd,frdf.pool,sum)[p.nms]
# gbly <- cbind(raw=rds.raw,flash=rds.flash,MID=rds.MID,primer=rds.dmult,
              # filt=rds.filt,ints=rds.ints)
# Tgbl <- colSums(gbly)
# pct.gbl <- round(Tgbl[-1]/Tgbl[-6]*100,2)

# pdf.flnm3 <- "GlobalYieldBarplots.pdf"
# pdf(file.path(repDir,pdf.flnm3),paper="a4",width=6.5,height=10)
# par(mfrow=c(2,1))

# pal1 <- brewer.pal(8,"Dark2")
# pal2 <- brewer.pal(8,"Pastel2")
# m <- ncol(gbly)
# ymx <- max(gbly)*1.2
# bp <- barplot(t(gbly),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              # ylim=c(0,ymx),xaxt="n",cex.axis=0.8,ylab="# reads")
# axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
# abline(h=0)	
# legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
# title(main="Yield on pools by analysis step")

# bp <- barplot(t(gbly/gbly[,1]*100),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              # ylim=c(0,117),xaxt="n",cex.axis=0.8,ylab="Percentage")
# axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
# grid(nx=NA,ny=NULL)	
# abline(h=0)	
# legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
# title(main="Yield on pools by analysis step")

# par(mar=c(5,8,4,6))

# bp <- barplot(pct.gbl,col="lavender",border="navy",ylim=c(0,max(pct.gbl)),
        # ylab="yield (%)")
# text(bp,50,pct.gbl,col="navy",font=2,cex=0.8)
# title(main="Global yield by step")		

# boxplot(frdf$Fn.rd,border="gray",ylab="# of reads",outline=FALSE,
        # ylim=range(frdf$Fn.rd))
# points(jitter(rep(1,nrow(frdf)),a=0.10),frdf$Fn.rd,pch="+",cex=0.8)
# title(main="Final coverage")

# dev.off()
			  
# txt.flnm2 <- "GlobalYield-SumRprt.txt"
# sink(file=file.path(repDir,txt.flnm2))

# cat("\n   Global yield by analysis step")
# cat("\n===================================\n")
# cat("\nIn number of reads:\n\n")
# gbly <- rbind(gbly,TOTAL=colSums(gbly))
# print(gbly)	
		
# cat("\nIn percentage by step:\n\n")
# print(round(gbly[,-1]/gbly[,-6]*100,2))

# cat("\nIn percentage referred to raw reads:\n\n")
# print(round(gbly/gbly[,1]*100,2))			
# sink()
