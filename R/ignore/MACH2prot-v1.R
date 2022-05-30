
###  Formatar sencers emplenant amb 0 per l'esquerra
######################################################
zeroFillInt2Char <- function(x,ln)
{ x <- paste("000000",x,sep="")
  substring(x,nchar(x)-ln+1)
}

###  Dóna nom als haplotipus en funció del nombre de mutacions respecte
###   la màster i de la seva frequencia poblacional, i salva les 
###   sequencies en un fitxer fasta.
########################################################################
SaveAaHaplotypes <- function(flnm,bseqs,nr)
{
  ##  Determinar dominant
  imstr <- which.max(nr)
  ##  Determinar diferències respecte la dominant
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[imstr])
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
  nms <- paste("Hpl",nm,zeroFillInt2Char(isq,4),sep=".")
  ##  Capçalera fasta amb nom, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")
  ##  Salvar a fasta
  write.fasta(AAStringSet(bseqs),flnm)
  list(bseqs=bseqs,nr=nr,nm=nm)
}

###  Recolapsar reads
######################
Recollapse <- function(seqs,nr)
{
  res <- tapply(nr,seqs,sum)
  nms <- names(res)
  res <- as.vector(res)
  names(res) <- NULL
  o <- order(res,decreasing=TRUE)
  list(bseqs=nms[o],nr=res[o])
}

###  Escriu unes seqüencies a un fitxer fasta
###############################################
write.fasta <- function(seqs,flnm)
{ writeXStringSet(seqs,flnm) }


###  Llegeix les seqüencies d'un fitxer fasta
###############################################
read.fasta <- function(flnm)
{ readAAStringSet(flnm) }

###  Aliniament múltiple per muscle de les seqs
###  Les seqüècies s'entren i es tornen com a DNAStringSet
############################################################
muscle <- "D:\\UltraSeq\\muscle\\muscle3.8.31_i86win32.exe"
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
          ignore.stdout=FALSE,invisible=FALSE)
  if( file.exists(res.file) )
    return( read.fasta(res.file) )
  return(NULL)
}

#---------------------------------------------------------------------------#

library(stringr)
library(Biostrings)
source("./R/seqanalfns.v4.5.R")

strip.gaps <- function(nts)
 paste(nts[nts!="-"],collapse="")

###  Fitxers a tractar i descripció
#####################################
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
FlTbl <- FlTbl[FlTbl$Str=="fw",]
o <- order(FlTbl$Pat.ID,FlTbl$Ampl.Nm,FlTbl$Str)
FlTbl <- FlTbl[o,]

flnms <- paste("MACHpl02",FlTbl$Pat.ID,FlTbl$Ampl.Nm,"fna",sep=".")
flnms <- file.path(mach.Dir,flnms)
pnms <- FlTbl$Pat.ID

p.flnms <- paste("MACH.aaP",FlTbl$Pat.ID,FlTbl$Ampl.Nm,"fna",sep=".")
p.flnms <- file.path(mach.P.Dir,p.flnms)
s.flnms <- paste("MACH.aaS",FlTbl$Pat.ID,FlTbl$Ampl.Nm,"fna",sep=".")
s.flnms <- file.path(mach.S.Dir,s.flnms)

for(i in 1:length(flnms))
{ 
  ###  Carregar fasta
  if(!file.exists(flnms[i])) next
  lst <- read.ampl.seqs(flnms[i],mnr=1)
  nr <- lst$IDs$nseqs
  ###  Eliminar gaps en l'aliniament
  seqs <- as.matrix(DNAStringSet(lst$seqs))
  seqs <- apply(seqs,1,function(nts)  paste(nts[nts!="-"],collapse=""))

  ###  Traduir en pauta P 
  p.seqs <- substr(seqs,1,(nchar(seqs)%/%3)*3)
  p.seqs <- as.character(translate(DNAStringSet(p.seqs)))
  lst2 <- Recollapse(p.seqs,nr)
  ###  Realinear
  p.seqs <- as.character(doMuscle(AAStringSet(p.seqs)))
  nr <- sapply(names(p.seqs),function(str) strsplit(str,"\\|")[[1]][2])
  nr <- as.integer(nr)
  ###  Salvar fasta
  lst3 <- SaveAaHaplotypes(p.flnms[i],p.seqs,nr)

  ###  Traduir en pauta S 
  s.seqs <- substring(seqs,2)
  s.seqs <- substr(s.seqs,1,(nchar(s.seqs)%/%3)*3)
  s.seqs <- as.character(translate(DNAStringSet(s.seqs)))
  lst2 <- Recollapse(s.seqs,nr)
  ###  Realinear
  s.seqs <- as.character(doMuscle(AAStringSet(s.seqs)))
  nr <- sapply(names(s.seqs),function(str) strsplit(str,"\\|")[[1]][2])
  nr <- as.integer(nr)
  ###  Salvar fasta
  lst3 <- SaveAaHaplotypes(s.flnms[i],s.seqs,nr)
}
