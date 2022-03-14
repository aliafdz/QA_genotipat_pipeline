# Descripció del fitxer QA inicial:
### Accept max of 5% bases below Q30 by read
ThrQ30 <- 0.05 # Aquesta variable s'hauria d'indicar en l'argument de la 
               # funció resultant d'aquest pas

library(ShortRead)
library(Biostrings)
library(stringr)

### Llegeix l'estructura de descripció de mostres
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      colClasses="character",stringsAsFactors=F)
### Llista de pools en samples
pools <- unique(samples$Pool.Nm)
### Llegeix descriptors de primers
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                      stringsAsFactors=F)

### Fitxers que resulten de Flash				
flnms <- list.files(flashDir)
# Guarda del fitxer només el nom del pool seguit de S1 o S2
snms <- sub("_flash\\.fastq$","",flnms)
# Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

# Genera dos fitxers amb extensió fastq (1 per pool) que resultaran d'aquest script
oflnms <- paste(parts[,"PatID"],"flashFilt.fastq",sep="_")
# Indica la ruta dels fitxers en la carpeta flashFilt
oflnms <- file.path(flashFiltDir,oflnms)
# Guarda la ruta dels fitxers inclosos en la carpeta Flash
flnms <- file.path(flashDir,flnms)

### Carrega el fitxer RData del rendiment de flash
load(file.path(repDir,"FLASH_table.RData"))

### Loop sobre pools
# Vector que inclou tants 0s com nº de pools
freads <- integer(length(snms))
for(i in 1:length(snms))
{ # Defineix les variables inicialment a 0
  raw.rds <- filt.rds <- 0
  ### Aplica streamer (iteració) en el fitxer fastq avaluat
  strm <- FastqStreamer(flnms[i],n=1e6)
  
  # Si existeix el fitxers flashFilt a la carpeta, l'elimina (ja que s'ha de generar ara)
  if(file.exists(oflnms[i]))
    file.remove(oflnms[i])
  appfl <- "w"
  ### Carrega el fitxer fastq per chuncks
  while(length(sqq <- yield(strm)))
  { ### Actualitza el nº de reads a cada iteració
    raw.rds <- raw.rds+length(sqq)
    ### Calcula els Phred scores (Codi ASCII-33) amb la funció 'quality()'
    phrsc <- as(quality(sqq),"matrix")
    ### A la matriu amb els Phred scores aplica el sumatori de les bases amb puntuació menor a 30 (per sota de Q30)   
    nl30 <- apply(phrsc,1,function(x) sum(x<30,na.rm=TRUE))
    ### Longituds de seqüència
    sqln <- width(sqq) # Cicles de seqüenciació	
    ### Divideix les bases<30 entre el nº de cicles -> fracció de les bases per read
  	fnl30 <- nl30/sqln
  	### Aplica el filtre per eliminar els reads amb més del 5% de les bases<Q30
  	sqq <- sqq[fnl30<=ThrQ30]
  	# Actualitza el nº de reads després de filtrar per Q30
  	filt.rds <- filt.rds+length(sqq)
  	
  	#writeXStringSet(seqs,oflnms[i],append=appfl) -> No sé si fa falta
  	# Genera el fitxer amb extensió fastq que es guardarà a flashFilt
  	writeFastq(sqq,oflnms[i],mode=appfl,compress=TRUE)
  	appfl <- "a"
  }
  close(strm)
  # Guarda per cada pool el nº de reads després del filtrat
  freads[i] <- filt.rds  
}
# Afegeix al fitxer RData que tenia de FLASH els valors de nº de reads després de filtrat Q30
flash.res$FiltQ30 <- freads
# També afegeix al fitxer el total de reads (resultat de sumar els que van fer extensió a FLASH i els que no)
flash.res$Raw <- flash.res$Extended+flash.res$NoExtd
# Actualitza el fitxer RData
save(flash.res,file=file.path(repDir,"FLASH_table.RData"))

# Genera una matriu amb les dades del fitxer RData: total de reads, extensió per FLASH i filtrats per Q30
res <- data.matrix(flash.res[,c("Raw","Extended","FiltQ30")])

# Genera el fitxer pdf on es generaran els gràfics
pdf(file.path(repDir,"FiltQ30.barplot.pdf"),paper="a4",width=6.6,height=10)
par(mfrow=c(2,1))
library(RColorBrewer)
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
# Defineix el límit de l'eix y a partir del valor màxim de la taula res
ymx <- max(res)*1.15
# Genera un gràfic de barres amb la trasposada de la taula res -> nº de reads per pas
bp <- barplot(t(res),col=pal2[1:3],border=pal1[1:3],beside=TRUE,
               xaxt="n",ylim=c(0,ymx))
axis(1,at=colMeans(bp),rownames(res),cex.axis=0.8,las=2)
legend("top",horiz=TRUE,fill=pal2,cex=0.8,legend=colnames(res))
# Calcula el % de reads respecte el total per veure el rendiment
resy <- round(res/res[,1]*100,1)
ymx <- 115
# Gràfic de barres que mostra el % de reads per pas de filtrat
bp <- barplot(t(resy),col=pal2[1:3],border=pal1[1:3],beside=TRUE,
               xaxt="n",ylim=c(0,ymx))
axis(1,at=colMeans(bp),rownames(res),cex.axis=0.8,las=2)
legend("top",horiz=TRUE,fill=pal2,cex=0.8,legend=colnames(res))
# Afegeix l'etiqueta dels % calculats al gràfic
y <- min(resy)/2
text(x=as.vector(bp),y=y,lab=as.vector(t(resy)),cex=0.5,font=2,srt=90,
     col=pal1[1:3])
dev.off()
# Guarda els resultats de la taula actualitzada en un fitxer .txt
sink(file.path(repDir,"FiltQ30_report.txt"))
print(flash.res)
sink()

