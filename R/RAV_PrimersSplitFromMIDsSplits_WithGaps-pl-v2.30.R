
###################################################################
###  AQUESTA VERSIÓ RETALLA ELS PRIMERS PER AMBDOS EXTREMS 
###   ANANT-LOS A BUSCAR, SENSE PREFIXAR L'AMPLADA DE L'AMPLICÓ.
###################################################################

library(Biostrings)
library(ShortRead)

#-----------------------------------------------------------------------#

###  Formatar sencers emplenant amb 0 per l'esquerra
######################################################
zeroFillInt2Char <- function(x,ln)
{ x <- paste("000000",x,sep="")
  substring(x,nchar(x)-ln+1)
}

###  Guarda unes seqüencies en format fasta
#############################################
write.fasta <- function(seqs,flnm)
{
  writeXStringSet(seqs,flnm)
}

###  Dóna nom als haplotipus en funció del nombre de mutacions respecte
###   la referència (que pot ser la màster) i de la seva frequencia
###   poblacional, i salva les sequencies en un fitxer fasta.
########################################################################
SaveAllHaplotypes <- function(bseqs,nr,flnm)
{
  code <- "Hpl"
  ##  Determinar diferències respecte la màster
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[1])
  nm <- nmismatch(psa)
  tnm <- table(nm)

  ##  Ordenar per nombre de mutacions
  o <- order(nm)
  bseqs <- bseqs[o]
  nr <- nr[o]
  nm <- nm[o]

  ##  Numero d'ordre dins de cada nombre de mutacions
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
  nms <- paste(code,nm,zeroFillInt2Char(isq,4),sep=".")

  ##  Capçalera fasta amb nom, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ##  Salvar a fasta
  write.fasta(DNAStringSet(bseqs),file.path(trimDir,flnm))

  list(bseqs=bseqs,nr=nr,nm=nm)
}

#-----------------------------------------------------------------------#

###  Llegim l'estructura de descripció de mostres
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      stringsAsFactors=F)
###  Llegim fitxers globals de HCV
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                    stringsAsFactors=F)

### Aquí les posicions de tall son estrictament les que deixen l'amplicó 
###   net de primers, i cal assegurar-ho!
primers$FW.tpos <- primers$FW.pos+nchar(primers$Primer.FW)
primers$RV.tpos <- primers$RV.pos-nchar(primers$Primer.RV)


###  Inicialitzacions
Ns <- nrow(samples)
Ns <- Ns*2					 
pr.res <- matrix(0,nrow=Ns,ncol=4)
colnames(pr.res) <-  c("Tot.reads","matches","shorts","fn.reads")
					 
FlTbl <- data.frame(File.Name=character(Ns),Pat.ID=character(Ns),
                    Ampl.Nm=character(Ns),Pr.ID=integer(Ns),
					Str=character(Ns),Pos=integer(Ns),Len=integer(Ns),
					Reads=integer(Ns),Hpls=integer(Ns),
					stringsAsFactors=FALSE)
					
###  Loop sobre pools en samples
k <- 0
pools <- unique(samples$Pool.Nm)
p.cv <- p.ok <- numeric(length(pools))
names(p.cv) <- names(p.ok) <- pools

sink(file.path(repDir,"AmpliconLengthsRprt.txt"))
pdf(file.path(repDir,"AmpliconLengthsPlot.pdf"),paper="a4",
              width="6",height=10.5)
par(mfrow=c(2,1))
			  
for(i in 1:length(pools))
{
  ###  Identificar mostres del pool
  idx <- which(samples$Pool.Nm==pools[i])
  ###  Fitxers de MIDs en pool
  flnms <- paste(pools[i],".MID",samples$MID[idx],".fna",sep="")
  
  ###  Loop sobre mostres en pool
  for(j in 1:length(idx))
  { jj <- idx[j]
  
    ###  Load fastq file
    if(! file.exists(file.path(splitDir,flnms[j]))) next 
	seqs <- readDNAStringSet(file.path(splitDir,flnms[j]))
	p.cv[i] <- p.cv[i] + length(seqs)
    ###  Get primer index
    ipr <- grep(samples$Primer.ID[jj],primers$Ampl.Nm)

	tt <- paste("Pool",pools[i],"  MID",samples$MID[jj],
	            "  Ampl",primers$Ampl.Nm[ipr])
	cat("\n\n",tt,"\n",sep="")
    		
    ###  Discard too short seqs
	seqs <- seqs[width(seqs)>180]

    ###  primer up matches 
    #################################
    pr.up <- primers$Primer.FW[ ipr ]
    up.matches <- vmatchPattern(pattern=pr.up,
                            subject=subseq(seqs,start=target.io,end=target.in),
                            max.mismatch=max.prdif,fixed=FALSE)
    delta <- target.io-1
    flags <- elementLengths(up.matches)>=1
    nr <- integer()
	shorts <- integer()
	trim.len <- 0
	
    up.flnm <- paste(samples$Patient.ID[jj],primers$Ampl.Nm[ipr],"PrFW.fna",
	                 sep=".")
	seqs.up <- ""
    if(sum(flags))
    { ###  Aquest pas és delicat, startIndex() té NULLs a la llista.
      pos <- sapply(startIndex(up.matches)[flags],function(x) x[[1]])+delta
      ###  Matches
      seqs.up <- seqs[flags]
	  ###  Update remaining seqs
	  seqs <- seqs[!flags]
	  ###  Trim 5' primer
      st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])   ###
      seqs.up <- subseq(seqs.up,start=st,end=width(seqs.up))   ###

      ###  Localitzar primer per l'altre extrem
      pr.p3 <- as.character(
	              reverseComplement(DNAString(primers$Primer.RV[ipr])))
      endp <- width(seqs.up)
      fstp <- endp-target.in
      delta <- fstp-1
      p3.matches <- vmatchPattern(pattern=pr.p3,
                            subject=subseq(seqs.up,start=fstp,end=endp),
                            max.mismatch=max.prdif,fixed=FALSE)
      flags <- elementLengths(p3.matches)>=1
      
      ###  Els que passen aquí tenen l'amplicó sencer
	  if(sum(flags))
	  { 
	    ###  Retallar-lo deixant l'amplicó net
	    seqs.up <- seqs.up[flags]
	    pos <- sapply(startIndex(p3.matches)[flags],function(x) x[[1]])+
	                  delta[flags]                             ###
		seqs.up <- subseq(seqs.up,start=1,end=pos-1)           ###
		
	    ###  Collapse sequences to haplotypes+frequencies
        sqtbl <- sort(table(as.character(seqs.up)),decreasing=TRUE)
        bseqs <- names(sqtbl)                
        names(bseqs) <- 1:length(bseqs)
        nr <- as.integer(sqtbl)
        lst <- SaveAllHaplotypes(bseqs,nr,up.flnm)
		nr <- lst$nr
		
		cat("\nForward seqs, table of read lengths (over 10 rd)\n")
        tbl.len <- tapply(nr,nchar(bseqs),sum)
		print(tbl.len[tbl.len>=10])
		plot(as.integer(names(tbl.len)),tbl.len,type="h",
		     xlab="Read length",ylab="# reads")
        title(main=paste(tt," Str FW"))			 
      }
    }
	k <- k+1
	FlTbl$File.Name[k] <- up.flnm
	FlTbl$Pat.ID[k] <- samples$Patient.ID[jj]
	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
	FlTbl$Pr.ID[k] <- ipr
	FlTbl$Str[k] <- "fw"
	FlTbl$Pos[k] <- primers$FW.tpos[ipr]
	FlTbl$Len[k] <- mean(width(seqs.up))
    FlTbl$Reads[k] <- sum(nr)				   
    FlTbl$Hpls[k] <- length(nr)				   

	p.ok[i] <- p.ok[i] + sum(nr)

	pr.res[k,1] <- length(seqs)
    pr.res[k,2] <- sum(flags)
    pr.res[k,3] <- sum(shorts)
    pr.res[k,4] <- sum(nr)
  
    ###  primer dn matches 
    #################################
    pr.dn <- primers$Primer.RV[ ipr ]
    dn.matches <- vmatchPattern(pattern=pr.dn,
                            subject=subseq(seqs,start=target.io,end=target.in),
                            max.mismatch=max.prdif,fixed=FALSE)
    delta <- target.io-1
    flags <- elementLengths(dn.matches)>=1
    nr <- integer()
    shorts <- integer()
	
    dn.flnm <- paste(samples$Patient.ID[jj],primers$Ampl.Nm[ipr],"PrRV.fna",
	                 sep=".")
	seqs.dn <- ""
    if(sum(flags))
    { ###  Matches
      seqs.dn <- seqs[flags]
	  ###  Update remaining seqs
	  seqs <- seqs[!flags]
	  ###  Trim primer
	  pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta
      st <- pos + (primers$RV.pos[ipr]-primers$RV.tpos[ipr])   ###
      seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))   ###
	  seqs.dn <- seqs.dn[width(seqs.dn)>target.in+5]
	  
      ###  Reverse complement down matches
      seqs.dn <- reverseComplement(seqs.dn)
      ###  ... and match up primer
      io <- max(target.io-5,1)
      delta <- io-1
      dn.matches <- vmatchPattern(pattern=pr.up,
                              subject=subseq(seqs.dn,start=io,end=target.in+5),
                              max.mismatch=max.prdif,fixed=FALSE)
      flags <- elementLengths(dn.matches)>=1

      ###  Els que passen aquí tenen l'amplicó sencer
      if(sum(flags))
	  { ###  Trim reversed dn seqs
        seqs.dn <- seqs.dn[flags]
        pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta
        st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])   ###
        seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))   ###
		
		###  Collapse sequences to haplotypes+frequencies
        sqtbl <- sort(table(as.character(seqs.dn)),decreasing=TRUE)
        bseqs <- names(sqtbl)                
        names(bseqs) <- 1:length(bseqs)
        nr <- as.integer(sqtbl)
	    lst <- SaveAllHaplotypes(bseqs,nr,dn.flnm)
        nr <- lst$nr

        cat("\nReverse seqs, table of read lengths (over 10 rd)\n")
        tbl.len <- tapply(nr,nchar(bseqs),sum)
        print(tbl.len[tbl.len>=10])
		plot(as.integer(names(tbl.len)),tbl.len,type="h",
		     xlab="Read length",ylab="# reads")
        title(main=paste(tt," Str RV"))			 
	  }
	}  

    k <- k+1
	FlTbl$File.Name[k] <- dn.flnm
	FlTbl$Pat.ID[k] <- samples$Patient.ID[jj]
	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
	FlTbl$Pr.ID[k] <- ipr
	FlTbl$Str[k] <- "rv"
	FlTbl$Pos[k] <- primers$FW.tpos[ipr]
	FlTbl$Len[k] <- trim.len
    FlTbl$Reads[k] <- sum(nr)				   
    FlTbl$Hpls[k] <- length(nr)				   
	
	p.ok[i] <- p.ok[i] + sum(nr)

	pr.res[k,1] <- length(seqs)
    pr.res[k,2] <- sum(flags)
    pr.res[k,3] <- sum(shorts)
    pr.res[k,4] <- sum(nr)
  }
}
sink()
dev.off()

###  Sincronitzar structures
fl <- FlTbl$Pr.ID>0
FlTbl <- FlTbl[fl,]
pr.res <- pr.res[fl,]
anms <- paste(FlTbl$Pat.ID,FlTbl$Ampl.Nm,sep=".")

rownames(samples) <- paste(samples$Patient.ID,samples$Primer.ID,sep=".")

###  Plot results
library(RColorBrewer)
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")

fw.idx <- which(FlTbl$Str=="fw")
rv.idx <- which(FlTbl$Str=="rv")
mprres <- data.frame(PatID=FlTbl$Pat.ID[fw.idx],
                     PrimerID=FlTbl$Ampl.Nm[fw.idx],
                     Treads=pr.res[fw.idx,1],
					 Shorts=pr.res[fw.idx,3]+pr.res[rv.idx,3],
                     FW.match=pr.res[fw.idx,4],
                     RV.match=pr.res[rv.idx,4],
					 Fn.reads=pr.res[fw.idx,4]+pr.res[rv.idx,4],
                     stringsAsFactors=FALSE)
mres <- mprres
T.reads <- apply(mres[,5:6],2, function(x)
              tapply(x,mres$PatID,sum))
if(length(unique(mres$PatID))==1)
{ x <- matrix(T.reads,nrow=1)
  rownames(x) <- mres$PatID[1]
  colnames(x) <- names(T.reads)
  T.reads <- x
}  
  
pdf.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=11)
par(mfrow=c(2,1),mar=c(7,4,4,2)+0.1)

ymx <- max(rowSums(T.reads))*1.2
barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

res.mat <- mres[,5:6]
ymx <- max(res.mat)*1.2
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
dev.off()

yield <- p.ok/p.cv*100
PoolTbl <- data.frame(MIDReads=p.cv,PrimerReads=p.ok,Pct=yield)

txt.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.txt",sep="_")
sink(file.path(repDir,txt.flnm))
cat("\nTable of reads identified by primer\n\n")
print(mprres)
cat("\n")
print(FlTbl[,-1])
cat("\nTotal reads identified by patient\n\n")
print(T.reads)
cat("\nYield by pool\n\n")
print(PoolTbl)
sink()

###  Save tables

# FlTbl <- FlTbl[FlTbl$Reads>0, ]  # Erase rows of null files
save(FlTbl,PoolTbl,file=file.path(repDir,"SplittedReadsFileTable.RData"))
sink(file.path(repDir,"SplittedReadsFileTable.txt"))
print(FlTbl)
sink()


###  Plot on A4 horizontally
pdf.flnm <- paste(proj.nm,"SplitByPrimersOnFlash-hz.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10,height=6)
par(mar=c(7,4,4,2)+0.1)

ymx <- max(res.mat)*1.2
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))

ymx <- max(rowSums(T.reads))*1.2
bp <- barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx),xaxt="n")
axis(side=1,at=bp,rownames(T.reads),cex.axis=0.6,las=2)
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

par(mfrow=c(1,2))

ymx <- max(c(res.mat[,1],res.mat[,2]))
rownames(primers) <- primers$Ampl.Nm
reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="fw"],"Region"]
boxplot(res.mat[,1]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,1],
       pch="+",cex=0.8)
title(main="Primers identified on forward reads")	   

reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="rv"],"Region"]
boxplot(res.mat[,2]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,2],
       pch="+",cex=0.8)
title(main="Primers identified on reverse reads")	   

dev.off()
