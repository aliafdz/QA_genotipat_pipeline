
library(ShortRead)
library(Biostrings)
library(stringr)

# Funció exactament igual que en el fitxer QCbyPos, però ara s'aplicarà sobre els resultats del filtrat Q30
fn.fastq <- function(flnm,ln) # ln= longitud d'amplicó
{ # Defineix una matriu de dimensions 5xln de 0s
  fnm.q <- matrix(0,nrow=5,ncol=ln)
  # Vector de 5 0s
  fnm.l <- numeric(5)
  # Variable numèrica
  nrds <- numeric()
  # Variable de nombre enter
  all.ln <- integer()
  
  ### Aplica streamer (iteració) en el fitxer fastq
  strm <- FastqStreamer(flnm,n=1e6) # és el valor de chunck.sz definit en el fitxer principal
  ### Carrega el fitxer fastq per chuncks
  nchk <- 0
  while(length(sqq <- yield(strm)))
  { nchk <- nchk+1 # Nº de cada iteració (chunck del fastq avaluat)
    nrds[nchk] <- length(sqq) # Guarda el nº de reads del chunck avaluat
	  
    ###  Phred scores. Codi ASCII-33
    # Funció 'quality()' retorna el valor de qualitat dels strings
    # Funció 'as()' permet fer coerció del resultat a matriu
    phrsc <- as(quality(sqq),"matrix")
	  
    # Guarda el valor mínim entre la variable ln i les columnes de la matriu amb les qualitats dels strings 
    nc <- min(ln,ncol(phrsc))
    # De les columnes 1 al mínim assignat abans, aplica els quantils del vector a les columnes
    # de la matriu amb les qualitats phrsc, i ho multiplica pel total de reads
    # Cada fila es un quantil i cada columna es un cicle de seqüenciació
    fnm.q[,1:nc] <- fnm.q[,1:nc] + 
	               apply(phrsc,2,quantile,p=c(0.05,0.25,0.5,0.75,0.95),
	                     na.rm=TRUE)[,1:nc] * nrds[nchk]
    
    ###  Longituds de seqüència
    sqln <- width(sqq) # Cicles de seqüenciació (longitud dels reads)
    all.ln <- c(all.ln,sqln) # Vector amb el nº de cicles com a nombres enters
    # Aplica els quantils per guardar-los en la llista de 5 columnes de 0 (1 columna per quantil)
    # i guarda els resultats multiplicats per la longitud de la seq
    fnm.l <- fnm.l + quantile(sqln,p=c(0.05,0.25,0.5,0.75,0.95),
	                     na.rm=TRUE) * nrds[nchk]
  }
  close(strm)
  # Retorna una llista amb 3 matrius:
  # 1- Fracció dels quantils de phred score entre total de reads
  # 2- Nº de reads de cada quantil entre total de reads
  # 3- Nº de cicles de seqüenciació
  return(list(fvnq=fnm.q/sum(nrds),fvnl=fnm.l/sum(nrds),all.ln=all.ln))
}  

# Exactament el mateix gràfic per resultats de FLASH que hi ha en el fitxer PoolQCbyPos:
# Gràfic de QC per posició i SW dels reads després de filtrar per FLASH
plot.flash <- function(fvnm)
  # N de columnes
{ nc <- ncol(fvnm)
  # Transforma en llista els resultats de quantils de phred score entre total de reads
  bxp.dt <- list(stats=fvnm,names=1:ncol(fvnm)) # ncol(fvnm) es podria canviar per nc
  # Dibuixa el gràfic i les línies de cada quantil
  plot(1:nc,fvnm[2,],type="l",col="darkgreen",ylim=c(0,40),xaxt="n",
       xlab="Position",ylab="Phred score")
  lines(1:nc,fvnm[3,],type="l",lwd=2,col="darkgreen")
  lines(1:nc,fvnm[4,],type="l",col="darkgreen")
  lines(1:nc,fvnm[1,],type="l",lty=3,col="darkgreen")
  lines(1:nc,fvnm[5,],type="l",lty=3,col="darkgreen")
  axis(1,at=seq(10,nc,10),lab=seq(10,nc,10),las=2,cex.axis=0.8)
  abline(h=c(20,23,30,33),lty=4,col="gray")
  legend("bottomleft",lty=c(3,1,1,1,3),lwd=c(1,1,2,1,1),cex=0.8,
         legend=rev(c(0.05,0.25,0.5,0.75,0.95)),title="quantiles")
  legend("bottom",lwd=2,col="darkgreen",legend="Flash reads",horiz=TRUE)		 
}	

#-------------------------------------------------------------------------#

runDir <- "./run"
flashDir <- "./flash"
repDir <- "./reports"
dataDir <- "./data"
flashFiltDir <- "./flashFilt"

###  Llegeix l'estructura de descripció de mostres
# La taula conté columnes: Patient.ID, MID, Primer.ID, Region, RefSeq.ID, Pool.Nm
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      colClasses="character",stringsAsFactors=F)
###  Llista de pools en samples
pools <- unique(samples$Pool.Nm)

###  Llegeix els descriptors de primers
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                      stringsAsFactors=F)

###  Longitud de l'amplicó major en cada pool
max.len.in.pool <- function(p)
{ # Guarda els items de la taula de mostres que corresponen al pool avaluat
  idx <- which(samples$Pool.Nm==p)
  # Si no hi ha cap item que correspongui retorna 0
  if(length(idx)==0) return(0)
  # Guarda l'item de l'amplicó de la taula primers que es troba en el ID del primer en taula de mostres
  idx <- which(primers$Ampl.Nm %in% samples$Primer.ID[idx])
  # Calcula el màxim de la resta entre la posició del primer reverse de l'amplicó avaluat i la posició del forward
  max(primers$RV.pos[idx]-primers$FW.pos[idx]+1)  
}
# Aplica la funció en els dos pools per guardar la longitud de l'amplicó major
# Com en aquest cas totes les mostres del mateix pool fan servir els mateixos primers, només hi ha un amplicó per pool
pln <- sapply(pools,max.len.in.pool)					  

###  Amb MID+M13
#pln <- pln + 2*(20+10)

###  Fitxers que resulten de Flash després de filtrar per Q30				
flnms <- list.files(flashFiltDir)
snms <- sub("_flashFilt\\.fastq$","",flnms)

###  Sincronitza la longitud dels amplicons dels pools amb el nom d'aquests
pln <- pln[snms]

###  Loop sobre pools
binsz <- 10
for(i in 1:length(snms))
{ # Genera el pdf indicat pel pool avaluat
  pdf.flnm <- paste("PoolFiltQCbyPos",snms[i],"pdf",sep=".")
  # Genera la ruta on es guardarà el pdf (reports)
  pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10.5,height=6.5)
  par(mfrow=c(2,1),mar=c(3,4,1.5,2)+0.1)
  
  # Aplica la primera funció sobre els fastq resultants de flashFilt
  lst <- fn.fastq(file.path(flashFiltDir,flnms[i]),pln[i])
  # Matriu amb fracció dels quantils de phred score entre total de reads
  fvnm <- lst$fvnq
  # Genera el gràfic de la segona funció definida al principi
  plot.flash(fvnm)
  
  # Nº de columnes dels resultats de flashFilt
  nc <- ncol(fvnm)
  # Aplica SW: per cada bloc de 10 columnes fa la mitjana del phred score de les files
  fvnm <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm[,iwin:(iwin+9),drop=FALSE]))
  # Gràfic SW dels resultats flashFilt
  plot.flash(fvnm)

  par(mfrow=c(1,2))
  # Gràfic de longitud de reads del flashFilt
  lnfrq <- table(lst$all.ln) # taula del nº de cicles
  x <- as.integer(names(lnfrq))
  y <- as.vector(lnfrq)
  plot(x,y,type="h",xlab="Read length",ylab="Frequency")
  title(main="Read lengths (Flash reads)")
  par(new=T)
  plot(x,cumsum(y/sum(y)),type="l",ylim=c(0,1),
       axes=F,xlab=NA,ylab=NA,col="blue")
  axis(side=4,col="blue",col.axis="blue",at=seq(0,1,0.1),las=2)
  grid()
  
  dev.off()
}
