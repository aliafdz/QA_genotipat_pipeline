
library(ShortRead)
library(Biostrings)
library(stringr)

# Funció que es cridarà després 
fn.fastq <- function(flnm,ln=301) # ln és 301 per defecte, però en realitat és la longitud de l'amplicó
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
    sqln <- width(sqq) # Cicles de seqüenciació 
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

# Gràfic de QC per posició i SW dels reads originals abans de flash
# Les línies dels quantils indiquen que com a màxim el % indicat dels reads tenen aquell Phred score
plot.R1R2 <- function(fvnm1,fvnm2,snm,SW=FALSE)
  # Columnes de la taula (fracció dels quantils de phred score entre total de reads)
{ nc <- ncol(fvnm1)
  # Dibuixa el gràfic amb el phred score per posició amb les línies dels quantils per R1 
  plot(1:nc,fvnm1[2,],type="l",col="navy",ylim=c(0,40),xaxt="n",
       xlab="",ylab="Phred score")
  lines(1:nc,fvnm1[3,],type="l",lwd=2,col="navy")
  lines(1:nc,fvnm1[4,],type="l",col="navy")
  lines(1:nc,fvnm1[1,],type="l",lty=3,col="navy")
  lines(1:nc,fvnm1[5,],type="l",lty=3,col="navy")
  axis(1,at=seq(10,nc,10),lab=seq(10,nc,10),las=2,cex.axis=0.8)
  # El mateix per R2
  nc <- ncol(fvnm2)
  lines(1:nc,fvnm2[2,],type="l",col="maroon")
  lines(1:nc,fvnm2[3,],type="l",lwd=2,col="maroon")
  lines(1:nc,fvnm2[4,],type="l",col="maroon")
  lines(1:nc,fvnm2[1,],type="l",lty=3,col="maroon")
  lines(1:nc,fvnm2[5,],type="l",lty=3,col="maroon")
  # Línies grises per alguns valors phred score que puguin ser d'utilitat
  abline(h=c(20,23,30,33),lty=4,col="gray")
  legend("bottomleft",lty=c(3,1,1,1,3),lwd=c(1,1,2,1,1),cex=0.8,
         legend=rev(c(0.05,0.25,0.5,0.75,0.95)),title="quantiles")
  legend("bottom",lwd=2,col=c("navy","maroon"),legend=c("R1","R2"),horiz=TRUE)
  # Títols en funció de si es SW o no
  if(SW)
    title(main=paste("Quality profile by SW (size 10,step 1):",snm))
  else
    title(main=paste("Quality profile by position:",snm))
}

# Gràfic de QC per posició i SW dels reads després de filtrar per FLASH
plot.flash <- function(fvnm)
  # Nº de columnes
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

###  Llegeix l'estructura de descripció de mostres
# La taula conté columnes: Patient.ID, MID, Primer.ID, Region, RefSeq.ID, Pool.Nm
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      colClasses="character",stringsAsFactors=F)
###  Llista de pools en samples
pools <- unique(samples$Pool.Nm)

###  Llegeix els descriptors de primers
# La taula conté les columnes: Ampl.Nm, Region, Primer.FW, Primer.RV, FW.pos, RV.pos, FW.tpos, RV.tpos, Aa.ipos, Aa.lpos
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

### S'afegeix a la longitud de l'amplicó la longitud del MID i M13
pln <- pln + 2*(20+10)

### Fitxers que resulten de Flash				
flnms <- list.files(flashDir)
# Guarda del fitxer només el nom del pool seguit de S1 o S2
snms <- sub("_flash\\.fastq$","",flnms)
# Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

### Fitxers R1 i R2 originals, de la carpeta run
R1.flnms <- paste(snms,"_L001_R1_001.fastq.gz",sep="")
R2.flnms <- paste(snms,"_L001_R2_001.fastq.gz",sep="")

### Sincronitza la longitud dels amplicons dels pools amb el nom d'aquests
pln <- pln[parts[,"PatID"]]

### Loop sobre pools
binsz <- 10
for(i in 1:length(snms))
{ # Genera el pdf PoolQCbyPos del pool avaluat
  pdf.flnm <- paste("PoolQCbyPos",parts[i,1],"pdf",sep=".")
  # Genera la ruta on es guardarà el pdf (reports)
  pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10.5,height=6.5)
  par(mfrow=c(2,1),mar=c(3,4,1.5,2)+0.1)
  
  # Aplica la funció del principi sobre els fitxers R1 i R2 del pool
  # Això és pel gràfic d'abans de flash, avalua la qualitat dels reads inicials
  lst1 <- fn.fastq(file.path(runDir,R1.flnms[i]))
  # Guarda la taula amb fracció dels quantils de phred score entre total de reads
  fvnm1 <- lst1$fvnq
  # Aplica la funció de nou sobre R2
  lst2 <- fn.fastq(file.path(runDir,R2.flnms[i]))
  fvnm2 <- lst2$fvnq
  
  # Genera el gràfic a partir de la segona funció definida
  plot.R1R2(fvnm1,fvnm2,snms[i])
  
  # Aplica la primera funció sobre els fastq resultants de flash
  lst <- fn.fastq(file.path(flashDir,flnms[i]),pln[i]) # pln per guardar com argument ln la longitud de l'amplicó!
  # Matriu amb fracció dels quantils de phred score entre total de reads
  fvnm <- lst$fvnq
  # Genera el gràfic de la 3a funció definida
  plot.flash(fvnm)
  
  # Nº de columnes dels resultats de R1
  nc <- ncol(fvnm1)
  # Aplica SW: per cada bloc de 10 columnes fa la mitjana del phred score de les files
  fvnm1 <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm1[,iwin:(iwin+9),drop=FALSE]))
  # El mateix per R2
  nc <- ncol(fvnm2)
  fvnm2 <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm2[,iwin:(iwin+9),drop=FALSE]))
  
  # Gràfic SW (sliding window, per regions) pels fitxers originals de run
  plot.R1R2(fvnm1,fvnm2,snms[i],TRUE)
  
  # Nº de columnes dels resultats de flash
  nc <- ncol(fvnm)
  fvnm <- sapply(1:(nc-10),function(iwin) 
             rowMeans(fvnm[,iwin:(iwin+9),drop=FALSE]))
  # Gràfic SW pels resultats de flash
  plot.flash(fvnm)
  
  # Gràfic de la longitud dels reads abans i després de flash
  par(mfrow=c(1,2))  
  stats=cbind(lst1$fvnl,lst2$fvnl,lst$fvnl)
  colnames(stats) <- c("R1","R2","Flash")
  bxp.dt <- list(stats=stats,names=colnames(stats))
  bxp(bxp.dt,pars=list(boxfilll="lavender"),border="navy",
      ylab="Read length",las=2,ylim=c(0,max(stats)))
  title(main="Read length distributions")
  
  # Gràfic de longitud de reads del flash
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
