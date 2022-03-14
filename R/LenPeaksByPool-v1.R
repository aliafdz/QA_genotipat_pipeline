
library(ShortRead)
library(Biostrings)
library(stringr)

# Aquesta funció té el mateix nom que en els dos fitxers previs, però aquest cop és més curta
# S'assembla més a la funció de PoolQCbyPos
fn.fastq <- function(flnm,ln=301)
{ # Defineix una matriu de dimensions 5xln de 0s
  fnm.q <- matrix(0,nrow=5,ncol=ln)
  # Vector de 5 0s
  fnm.l <- numeric(5)
  # Variable numèrica
  nrds <- numeric()
  # Variable de nombre enter
  all.ln <- integer()
  
  ### Aplica streamer (iteració) en el fitxer fastq
  strm <- FastqStreamer(flnm,n=chunck.sz) # chunck.sz definit en el fitxer principal
  ### Carrega el fitxer fastq per chuncks
  while(length(sqq <- yield(strm)))
  { ###  Longituds de seqüència 
    all.ln <- c(all.ln,width(sqq)) # Vector amb el nº de cicles com a nombres enters
  }
  close(strm)
  # Retorna el vector
  return(all.ln)
}  

#-------------------------------------------------------------------------#

runDir <- "./run"
flashDir <- "./flash"
repDir <- "./reports"
dataDir <- "./data"

### Llegeix l'estructura de descripció de mostres
# La taula conté columnes: Patient.ID, MID, Primer.ID, Region, RefSeq.ID, Pool.Nm
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      colClasses="character",stringsAsFactors=F)
###  Llista de pools en samples
pools <- unique(samples$Pool.Nm)

###  Fitxers que resulten de Flash				
flnms <- list.files(flashDir)
# Guarda del fitxer només el nom del pool seguit de S1 o S2
snms <- sub("_flash\\.fastq$","",flnms)
# Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

# Genera el fitxer pdf indicat en el directori de reports
pdf(file.path(repDir,"PoolReadLengths.pdf"),paper="a4",width=6,height=10)
par(mfrow=c(2,1),mar=c(5,4,4,2.5)+0.1)
par(mfrow=c(2,1))

###  Loop sobre pools
for(i in 1:length(snms))
{ # Aplica la funció definida (FastqStreamer) en els fitxers de flash
  vln <- fn.fastq(file.path(flashDir,flnms[i]))
  # Taula per indicar la freqüència de cada nº de cicles extrets en vln
  lnfrq <- table(vln)
  # Defineix els valors del nº de cicles (longitud dels reads)
  x <- as.integer(names(lnfrq))
  # Defineix les freqüències
  y <- as.vector(lnfrq)
  # Gràfic de la longitud dels reads i la seva freq
  plot(x,y,type="h",xlab="Read length",ylab="Frequency")
  title(main=paste(snms[i],"- read lengths"),line=2)
  par(new=T)
  # Calcula la freq acumulada de la longitud dels reads per dibuixar la línea al gràfic
  plot(x,cumsum(y/sum(y)),type="l",ylim=c(0,1),
       axes=F,xlab=NA,ylab=NA,col="blue")
  # L'eix y de la línia de freq és de color blau amb intervals de 0.1
  axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2,cex.axis=0.8)
  grid()
  # Guarda en una variable les longituds de read amb freq acumulada major al 5%
  idx <- which( y/sum(y) >= 0.05 )
  x[idx]
  # Mostra un altre títol de gràfic indicant on es troben els pics de longitud de seq
  title(main=paste("Peaks at",paste(x[idx],collapse=", ")),
        cex.main=0.8,line=0.8)
}

dev.off()
