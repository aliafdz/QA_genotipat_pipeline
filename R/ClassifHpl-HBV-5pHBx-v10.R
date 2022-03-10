##############################################################################
###                    MODUL: ClassifHpl-HBV          
###
###  Classifica tots els haplotips de consens per genotips
###  Es genera un informe sobre els resultats obtesos.
##############################################################################

###  Guarda unes seqüencies en format fasta
#############################################
write.fasta <- function(seqs,flnm)
{
  writeXStringSet(DNAStringSet(seqs),flnm)
}


###  Aliniament múltiple per muscle de les seqs de referència
###    i les d'un pacient. 
###  Les seqüècies s'entren com a DNAStringSet, però es tornen
###    aliniades com a una matriu en format binari de 'ape' 
######################################################################
#muscle <- "D:/UltraSeq/muscle/muscle3.8.31_i86win32.exe"
#muscle <- "E:/UltraSeq/muscle/muscle3.8.31_i86win32.exe"
#muscle <- "C:/Docs/UltraSeq/muscle/muscle3.8.31_i86win32.exe"
muscle.cl.opts <- c("-log muscle.log")
doMuscle <- function(seqs,refs)
{ tmp.file <- file.path(tempDir,"muscleInFile.fna")
  res.file <- file.path(tempDir,"muscleOutFile.fna")
  if(file.exists(res.file)) file.remove(res.file)
  write.fasta( c(refs,seqs), tmp.file )
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
plot.Px.HCdend <- function(rseqs,tit,model="K80",gamma=FALSE,nchtp=3)
{ library(RColorBrewer)
  cls <- c(brewer.pal(8,"Dark2"),"black","red","blue")
  #cls <- c(brewer.pal(6,"Dark2"),"black")
  dm <- dist.dna(rseqs,model,gamma=gamma,as.matrix=FALSE)

  nms <- attr(dm,"Labels")
  tp <- as.factor(substr(nms,1,nchtp))

  hc <- hclust(as.dist(dm),method="average")
  dend_colo1 <- as.dendrogram(hc)
  for(i in 1:nlevels(tp))
   dend_colo1 <- dendrapply(dend_colo1,edgeCol, fgr=cls[i], 
                  lwd=2,keys=nms[tp==levels(tp)[i]]) 

  omar <- par(mar=c(3,2,2.5,5))
  par(cex.main=1,cex.lab=0.8,cex.axis=0.8)
  plot(dend_colo1,nodePar=list(pch=NA,lab.cex=0.6),horiz=TRUE)
  tt1 <- paste("UPGMA tree (",model,sep="",collapse="")
  if(gamma) tt1 <- paste(tt1,"; gamma=",gamma,sep="",collapse="")
  title(paste(tt1,"):  ",tit,sep="",collapse=""),cex.main=0.8,line=1)
  par(mar=omar)
  invisible(labels(dend_colo1))
}


###  Plot a MDS map, PC1/PC2, PC1/PC3 and PC2/PC3 RefSeqs 
###########################################################
plot.Px.mds <- function(rseqs,tit="Px",model="K80",
                          gamma=FALSE,nchtp=1,nmln=2)
{ library(RColorBrewer)
  cls <- c(brewer.pal(8,"Dark2"),"black","red","blue")
  #cls <- c(brewer.pal(6,"Dark2"),"black")

  rownames(rseqs) <- substr(rownames(rseqs),1,nmln)
  tp <- as.factor(substr(rownames(rseqs),1,nchtp))
  dm <- dist.dna(rseqs,model,gamma=gamma,as.matrix=FALSE)
  mds <- pcoa(dm, correction="lingoes", 
              rn=substr(rownames(data),1,5))

  tt1 <- paste("MDS map (",model,sep="",collapse="")
  if(gamma) tt1 <- paste(tt1,"; gamma=",gamma,sep="",collapse="")

  plot(mds$vectors[,1],mds$vectors[,2],type="n",asp=1,
       xlab="PC1",ylab="PC2")
  text(mds$vectors[,1],mds$vectors[,2],rownames(mds$vectors),cex=0.7,
     col=cls[as.integer(tp)],font=2)
  abline(h=0,v=0,lty=4,col="lightgray")
  title(paste(tt1,"):  ",tit,sep="",collapse=""),cex.main=1,line=1)

  plot(mds$vectors[,1],mds$vectors[,3],type="n",asp=1,
       xlab="PC1",ylab="PC3")
  text(mds$vectors[,1],mds$vectors[,3],rownames(mds$vectors),cex=0.7,
     col=cls[as.integer(tp)],font=2)
  abline(h=0,v=0,lty=4,col="lightgray")

  plot(mds$vectors[,2],mds$vectors[,3],type="n",asp=1,
       xlab="PC2",ylab="PC3")
  text(mds$vectors[,2],mds$vectors[,3],rownames(mds$vectors),cex=0.7,
     col=cls[as.integer(tp)],font=2)
  abline(h=0,v=0,lty=4,col="lightgray")
}


###  Sort haplotype names
order.hpl <- function(hnms)
{ nms <- sapply(hnms,function(hnm) strsplit(hnm,split="\\|")[[1]][1])
  mto <- t(sapply(nms,function(nm) 
               as.integer(strsplit(nm,split="\\.")[[1]][2:3])))
  order(mto[,1],mto[,2])
}

sort.hpl <- function(hnms)
{ 
  hnms[order.hpl(hnms)]
}

############################################################################
###  MAIN LINE  ###

library(stringr)
library(Biostrings)
library(ape)

hplDir  <- mach.Dir

###  Codi de la DB rule de Cuadras
source("./R/DBrule_Genotipat.v2.R")

## Noms d'amplico
ampl.nm <- "1255:1611"
a.nms <- ampl.nm

# MACHpl02.CCG0057.HBV.1234.1631.fna
flnms <- dir(hplDir,"MACHpl02.+\\.HBV\\.1234\\.1631\\.fna$")  
pnms <- sub("MACHpl02\\.","",flnms)
pnms <- str_extract(pnms,"^[A-Za-z0-9\\-]+")
snms <- rep(ampl.nm,length(pnms))
flnms <- file.path(hplDir,flnms)
n <- length(flnms)

###  Carregar sequencies de referencia
gen.nm <- tit <- c("HBV_RefSeqs_1255_1611")
rs.flnm <- paste(gen.nm,".fna",sep="")
rs.flnm <- file.path(dataDir,rs.flnm)
rseq <- as.character(readDNAStringSet(rs.flnm))
gnms <- sort(unique(substr(names(rseq),1,1)))

#--------------------------------------------------------------------------#

pdf.flnm <- file.path(resultsDir,"HBV_1255_1611_HplGenotyping.pdf")
pdf(file=pdf.flnm,paper="a4",width=5,height=11)

res <- data.frame(Pnm=character(),Ampl=character(),Hpl=character(),
                  Type=character(),MinPhi2=numeric(),NextMin=numeric(),
				  Ratio=numeric(),stringsAsFactors=FALSE)
res.dst <- data.frame(Phi2.1=numeric(),Phi2.2=numeric(),
                  Phi2.3=numeric(),Phi2.4=numeric(),Phi2.5=numeric(),
                  Phi2.6=numeric(),Phi2.7=numeric(),Phi2.8=numeric(),
                  Phi2.9=numeric(),Phi2.10=numeric(),				  
                  stringsAsFactors=FALSE)
res.dst <- res.dst[,1:length(gnms)]
colnames(res.dst) <- Phi.gnms <- paste("Phi2",gnms,sep=".")

k <- 0
for(i in 1:n)
{ 
  ###  Carregar haplotips del pacient
  flnm <- flnms[i]
  if(!file.exists(flnms[i])) next
  pseq <- as.character(readDNAStringSet(flnm))
  names(pseq) <- paste("x",names(pseq),sep="")

  ###  Eliminar Hpl < 1%
  nr <- sapply(names(pseq),function(str) 
                          as.integer(strsplit(str,split="\\|")[[1]][2])) 
  fl <- nr/sum(nr)*100 >= 1
  pseq <- pseq[fl]
  nr <- nr[fl]
  ###  ... i no màster per sota de 100 reads
  fl <- nr >= 100
  fl[1] <- TRUE
  pseq <- pseq[fl]
  
  ###  Aliniament múltiple de RefSeqs i haplotips del pacient
  bseqs <- doMuscle(pseq,rseq)
  onms <- rownames(bseqs)
  #rownames(bseqs) <- substr(onms,1,20)
  if( is.null(bseqs) ) next
 
  ###  Arbre filogenertic i MDS sobre els tres primers CP
  tt <- paste(pnms[i],snms[i])
  plot.Px.HCdend(bseqs,tt,model="K80",nchtp=1)
  opar=par(mfrow=c(3,1),mar=c(4,4,3,2)+0.1)
  plot.Px.mds(bseqs,tit=tt,model="K80",gamma=FALSE,nchtp=1,nmln=2)
  par(opar)

  ###  Classificació per la DB rule de Cuadras
  dm <- dist.dna(bseqs,model="K80",as.matrix=TRUE,pairwise.deletion=TRUE)
  jdx <- which(substr(rownames(bseqs),1,4)=="xHpl")
  dgrp <- dm[-jdx,-jdx]
  grp <- factor(substr(rownames(dgrp),1,1))
  hr <- as.integer(grp)
  k1 <- k+1
  for(ii in jdx)
  { d <- dm[ii,-jdx]
    dsc <- DBrule(dgrp, hr, d, my.names=levels(grp))
    k <- k+1
    res[k,1] <- pnms[i]
    res[k,2] <- snms[i]
    res[k,3] <- rownames(bseqs)[ii]
    res[k,4] <- dsc$Type[1]
    res[k,5] <- signif(min(dsc$Phi2),3)
    res[k,6] <- signif(sort(dsc$Phi2)[2],3)
    res[k,7] <- round(res[k,6]/res[k,5],3)
	
    res.dst[k,Phi.gnms] <- dsc$Phi2[Phi.gnms]
  }

  ###  Convertir a vector de strings
  sqs <- apply(as.character(bseqs[jdx,]),1,function(vch) paste(vch,collapse=""))

  ###  Imposar ordre
  o <- order.hpl(names(sqs))
  sqs <- sqs[o]
  hdnm <- paste((res[k1:k,4])[o],pnms[i],sep=".")
  names(sqs) <- paste(hdnm,substring(sort.hpl(onms[jdx]),6),sep=".")
  res[k1:k,] <- res[(k1:k)[o],]
  res.dst[k1:k,] <- res.dst[(k1:k)[o],]
}
dev.off()

res$Hpl <- sub("xHpl.","",res$Hpl)
res$Ratio <- round(res$Ratio,2)
res$flag <- ifelse(abs(res$Ratio)<2,"*","")

#--------------------------------------------------------------------------#

###  Expressio de resultats
sink(file.path(resultsDir,"HBV_1255_1611_HplGenotypingRprt.txt"))
cat("\n   GENOTYPING HBV CONSENSUS HAPLOTYPES")
cat("\n=========================================\n")
cat("\n Results of discriminant analysis by the DB rule:\n\n")
print(res)
cat("\n")
print(signif(res.dst,3))

###  Cobertures per haplotip i pacient
nr <- sapply(res$Hpl,function(str) 
                          as.integer(strsplit(str,split="\\|")[[1]][2])) 
covs <- tapply(nr,res$Pnm,sum)

###  Resultat global per pacient
gtps <- gnms
gtbl <- matrix(0,nrow=length(covs),ncol=length(gtps))
colnames(gtbl) <- gtps
rownames(gtbl) <- names(covs)

for(i in 1:length(covs))
{ pid <- names(covs)[i]
  fl <-  res$Pnm==pid
  tbl1 <- tapply(nr[fl],res$Type[fl],sum)
  gtbl[i,names(tbl1)] <- tbl1
}

gtbl2 <- t(apply(gtbl,1,function(x) round(x/sum(x)*100,2)))

df.gtbl <- data.frame(gtbl)
df.gtbl$flag <- ""
df.gtbl2 <- data.frame(gtbl2)
df.gtbl2$flag <- ""
for(i in 1:nrow(df.gtbl))
{ if(any(res[res$Pnm==rownames(df.gtbl)[i],"flag"]=="*"))
  { df.gtbl$flag[i] <- "*"
    df.gtbl2$flag[i] <- "*"
  }
}

cat("\n\nCoverage by patient:\n")
print(covs)
cat("\n\nSummary by coverage in number of reads classified:\n")
print(df.gtbl)
cat("\n\nSummary in percentage by genotype:\n")
print(df.gtbl2)

sink()

#--------------------------------------------------------------------------#

sink(file.path(resultsDir,"HBV_1255_1611_HplGenotyping_SummaryRprt.txt"))
cat("\n   GENOTYPING HBV CONSENSUS HAPLOTYPES")
cat("\n=========================================\n")
cat("\n Results of discriminant analysis by the DB rule:\n\n")
cat("\n\nCoverage by patient:\n")
print(covs)
cat("\n\nSummary by coverage in number of reads classified:\n\n")
print(df.gtbl)
cat("\n\nSummary in percentage by genotype:\n\n")
print(df.gtbl2)
cat("\n\nGlobally by genotype:\n\n")
print(round(apply(gtbl2,2,sum)/nrow(gtbl2),2))
cat("\nsorted in decreasing order:\n\n")
print(sort(round(apply(gtbl2,2,sum)/nrow(gtbl2),2),decreasing=TRUE))

sink()

#--------------------------------------------------------------------------#

library(tidyverse)

P <- gtbl2 %>% data.frame() %>%
  rownames_to_column(var='ID') %>%
  pivot_longer(-ID,names_to='Genotype',values_to='Freq') %>%
  ggplot() +
    geom_col(aes(x=ID,y=Freq,fill=Genotype),col='lavender') +
    scale_fill_brewer(palette = "Dark2") +
    labs(x='',y='Frequency (%)',title='Genotyping HDV') +
    theme_bw(base_size=14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="bottom")

pdf.flnm <- file.path(resultsDir,'HBV_1255_1611_GenotypingBarplot.pdf')
pdf(pdf.flnm,paper='a4',width=7,height=6)
print(P)
dev.off()


P <- res %>%
  select(Pnm,Type,Ratio) %>%
  mutate(Ratio=ifelse(Ratio > 10 | Ratio < 0,10,Ratio)) %>%
  ggplot() +
    geom_boxplot(aes(x=Pnm,y=Ratio),col='gray',alpha=0.3,outlier.size=0)+
    geom_jitter(aes(x=Pnm,y=Ratio,col=Type),width=0.2,height=0,size=0.7) +
    geom_hline(aes(yintercept=2),linetype=4,col='maroon')+
    labs(x='',y=expression(paste("Distance ratio  " * bold(phi[2]^2/phi[1]^2))),
         title='Haplotype genotyping quality') +
    ylim(0,NA) +
    theme_bw(base_size=14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position="bottom")
  
pdf.flnm <- file.path(resultsDir,'HBV_1255_1611_HplGenotypingConfidence.pdf')
pdf(pdf.flnm,paper='a4',width=7,height=6)
print(P)
dev.off()

#--------------------------------------------------------------------------#

save(covs,res,res.dst,gtbl,gtbl2,df.gtbl,df.gtbl2,
     file=file.path(resultsDir,"5pHBx_Classif.RData"))

#--------------------------------------------------------------------------#
