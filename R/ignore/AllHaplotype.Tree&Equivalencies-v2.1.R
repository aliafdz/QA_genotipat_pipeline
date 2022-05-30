##############################################################
###   GENOTYPING VIRUS SEQUENCES FROM SANGER SEQUENCES     ###
##############################################################

library(stringr)
library(Biostrings)
library(ape)
source("./R/seqanalfns.v4.5.R")

###  Guarda unes seqüencies en format fasta
#############################################
write.fasta <- function(seqs,flnm)
{
  writeXStringSet(seqs,flnm)
}

###  Llegeix seqüencies en format fasta
#############################################
read.fasta <- function(flnm)
{
  readDNAStringSet(flnm)
}


###  Aliniament múltiple per muscle de les seqs  
###  Les seqüècies s'entren com a DNAStringSet, però es tornen
###    aliniades com a una matriu en format binari de 'ape' 
######################################################################
#muscle <- "C:/Muscle/muscle3.8.31_i86win32.exe"
#muscle <- "D:/UltraSeq/muscle/muscle3.8.31_i86win32.exe"
#muscle <- "C:/Docs/UltraSeq/muscle/muscle3.8.31_i86win32.exe"
muscle.cl.opts <- c("-log muscle.log")
doMuscle <- function(seqs)
{ tmp.file <- file.path(tempDir,"muscleInFile.fna")
  res.file <- file.path(tempDir,"muscleOutFile.fna")
  if(file.exists(res.file)) file.remove(res.file)
  write.fasta( seqs, tmp.file )
  in.file <- paste("-in ",tmp.file,sep="")
  out.file <- paste("-out ",res.file,sep="")
  command <- paste(muscle,in.file,out.file,
                paste(muscle.cl.opts,collapse=" "),sep=" ")
  system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
          ignore.stdout=FALSE,invisible=TRUE)
  if( file.exists(res.file) )
    return( read.dna(res.file,format="fasta") )
  return(NULL)
}



###  Highlight the edges bringing to a set of labels
######################################################
edgeCol <- function(dend, keys, fgr="red", lwd=1, ...) 
{ myattr <- attributes(dend)
  if(is.leaf(dend)) 
  { # Si es una fulla amb etiqueta en keys
    if(length(which(keys==myattr$label))==1)
      attr(dend,"edgePar") <- c(myattr$edgePar, 
           list(col=fgr,lwd=lwd))
  } else {
    lbl <- labels(dend)
    # Si es una branca amb totes les fulles en keys
    if( length(intersect(keys,lbl))==length(lbl) )
      attr(dend,"edgePar") <- c(myattr$edgePar,
           list(col=fgr,lwd=lwd))
  }
  return(dend)
}


###  Plot the HC dendrogram
#############################
plot.Px.HCdend <- function(data,tit,model="K80",gamma=FALSE,nch=2)
{ 
  library(RColorBrewer)

  # cls <- c(brewer.pal(8,"Dark2"),"black")
  cls <- rep("black",25)
  
  #rownames(data) <- substr(rownames(data),1,15)
  dm <- dist.dna(data,model,gamma=gamma,as.matrix=FALSE,pairwise.deletion=TRUE)

  #nv28 <- substr(attr(dm,"Labels"),1,15)
  #attr(dm,"Labels") <- nv28
  #tp <- as.factor(substr(nv28,1,nch))
  nv28 <- rownames(data)
  tp <- as.factor(substr(nv28,1,nch))

  hc <- hclust(as.dist(dm),method="average")
  dend_colo1 <- as.dendrogram(hc)
  for(i in 1:nlevels(tp))
   dend_colo1 <- dendrapply(dend_colo1,edgeCol, fgr=cls[i], 
                  lwd=2,keys=nv28[tp==levels(tp)[i]]) 

  omar <- par(mar=c(3,2,2.5,6))
  par(cex.main=1,cex.lab=0.8,cex.axis=0.8)
  plot(dend_colo1,nodePar=list(pch=NA,lab.cex=0.5),horiz=TRUE)
  abline(v=seq(0,80,10),lty=4,col="lightgray")
  gtnm <- levels(tp)
  #legend("topleft",legend=gtnm,lwd=2,cex=0.7,bg="white",
  #       col=cls[1:nlevels(tp)]) 
  tt1 <- paste("UPGMA tree (",model,sep="",collapse="")
  if(gamma) tt1 <- paste(tt1,"; gamma=",gamma,sep="",collapse="")
  title(paste(tt1,"):  ",tit,sep="",collapse=""),cex.main=1,line=1)
  par(mar=omar)	 
}


#####################
###   MAIN LINE   ###
#####################

###  Models d'evolució genètica per les ref.seqs
##################################################
model <- "N"
gamma <- FALSE
###  Nombre de caracters en nov28 haplotips pel colorejat de branques
nch <- 1

ntDir <- mach.Dir

#--------------------------------------------------------------------------#

tt <- "BQ43_HBV.1234.1631"
flnms <- list.files(path=ntDir,patt="^MACHpl02\\..*\\.HBV.1234.1631.fna$")
snms <- sub("MACHpl02\\.","",flnms)
snms <- str_extract(snms,"^[0-9A-Za-z_\\-]+")

if(length(snms) > 1) ### Cal que hi hagi més d'un pacient.
{ ###  Agregar sequencies de totes les mostres en l'experiment
  all.seqs <- DNAStringSet()
  for(i in 1:length(snms))
  { seqs <- readDNAStringSet(file.path(ntDir,flnms[i]))
    names(seqs) <- sub("Hpl",snms[i],names(seqs))
    all.seqs <- c(all.seqs,seqs)
  }

  ###  Aliniament múltiple amb muscle
  aseqs <- doMuscle(all.seqs)
  fna.flnm <- paste(tt,"AllSamples.fna",sep=".")
  file.copy("muscleOutFile.fna",file.path(resultsDir,fna.flnm),overwrite=TRUE)

  ###  Arbre UPGMA amb tots els haplotips  
  pdf.flnm <- paste(tt,"AllSamples.HplUPGMATree.pdf",sep=".")
  pdf(file=file.path(resultsDir,pdf.flnm),width=6,height=length(all.seqs)/5)
  plot.Px.HCdend(aseqs,tt,model=model,gamma=gamma,nch)
  dev.off()

  ###  Equivalències (contaminacions ?)
  ###  Cal que siguin estrictament iguals per registrar-les 
  fl <- rep(TRUE,nrow(aseqs))
  k <- 1
  eq.lst <- list()
  while(sum(fl)>1)
  { idx <- which(fl)
    i <- idx[1]
    fl[i] <- FALSE
    idx <- idx[-1]
    fl2 <- sapply(idx,function(j) all(aseqs[i,]==aseqs[j,]))
    if(sum(fl2))
    { eq.lst[[k]] <- rownames(aseqs)[c(i,idx[fl2])]
      k <- k+1
	  fl[idx[fl2]] <- FALSE
    }
  }

  txt.flnm <- paste(tt,"HplEquivalencies.txt",sep=".")
  sink(file.path(resultsDir,txt.flnm))
  print(eq.lst[order(sapply(eq.lst,length),decreasing=TRUE)])
  sink()
}

#--------------------------------------------------------------------------#

tt <- "BQ43_HBV.2822.3296"
flnms <- list.files(path=ntDir,patt="^MACHpl02\\..*\\.HBV.2822.3296.fna$")
snms <- sub("MACHpl02\\.","",flnms)
snms <- str_extract(snms,"^[0-9A-Za-z_-]+")

if(length(snms) > 1) ### Cal que hi hagi més d'un pacient.
{ ###  Agregar sequencies de totes les mostres en l'experiment
  all.seqs <- DNAStringSet()
  for(i in 1:length(snms))
  { seqs <- readDNAStringSet(file.path(ntDir,flnms[i]))
    names(seqs) <- sub("Hpl",snms[i],names(seqs))
    all.seqs <- c(all.seqs,seqs)
  }

  ###  Aliniament múltiple amb muscle
  aseqs <- doMuscle(all.seqs)
  fna.flnm <- paste(tt,"AllSamples.fna",sep=".")
  file.copy("muscleOutFile.fna",file.path(resultsDir,fna.flnm),overwrite=TRUE)

  ###  Arbre UPGMA amb tots els haplotips  
  pdf.flnm <- paste(tt,"AllSamples.HplUPGMATree.pdf",sep=".")
  pdf(file=file.path(resultsDir,pdf.flnm),width=6,height=length(all.seqs)/5)
  plot.Px.HCdend(aseqs,tt,model=model,gamma=gamma,nch)
  dev.off()

  ###  Equivalències (contaminacions ?)
  ###  Cal que siguin estrictament iguals per registrar-les 
  fl <- rep(TRUE,nrow(aseqs))
  k <- 1
  eq.lst <- list()
  while(sum(fl)>1)
  { idx <- which(fl)
    i <- idx[1]
    fl[i] <- FALSE
    idx <- idx[-1]
    fl2 <- sapply(idx,function(j) all(aseqs[i,]==aseqs[j,]))
    if(sum(fl2))
    { eq.lst[[k]] <- rownames(aseqs)[c(i,idx[fl2])]
      k <- k+1
	  fl[idx[fl2]] <- FALSE
    }
  }

  txt.flnm <- paste(tt,"HplEquivalencies.txt",sep=".")
  sink(file.path(resultsDir,txt.flnm))
  print(eq.lst[order(sapply(eq.lst,length),decreasing=TRUE)])
  sink()
}

#--------------------------------------------------------------------------#
