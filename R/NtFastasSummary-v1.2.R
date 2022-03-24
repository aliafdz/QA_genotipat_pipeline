
library(Biostrings)
library(stringr)

### Funció per llegir les seqüències aliniades dels amplicons
  # Els arguments corresponen al nom del fitxer on es troben les seqs i el valor
  # de mínim nº de reads assignat (per defecte 1)
read.ampl.seqs <- function(flnm,mnr=1) # flnm= file.path(mach.Dir,flnms[i])
{ # Llegeix el fitxer fasta de l'argument i el transforma en un objecte DNAStringSet
  seqs <- as.character(readDNAStringSet(flnm))
  # Guarda els noms de les seqüències del fitxer, és a dir dels haplotips
  # Aquests haplotips són els que han fet intersecció entre ambdues cadenes
  IDstr <- names(seqs)
  # Guarda el total de seqs del fitxer (nº d'haplotips coincidents)
  n <- length(IDstr)
  
  # Variable de caràcter buida amb tantes entrades com haplotips trobi al fitxer 
  nms <- character(length=n)
  # Variable numèrica buida amb tantes entrades com haplotips trobi al fitxer
  nr <- integer(n)
  # Bucle sobre els haplotips del fitxer (coincidents)
  for(j in 1:n)
  { # Separa el nom de l'haplotip avaluat segons el símbol "|"
    strs <- strsplit(IDstr[j],split="\\|")[[1]]
    # Guarda a la variable buida d'abans el primer element del nom separat,
    # que correspon a la consecució de "Hpl", el nº de mutacions amb la seq màster
    # i el nº d'ordre dins del grup amb les mateixes mutacions
    nms[j] <- strs[1]
    # Si la variable strs (separació del nom) conté més d'un element,
    if(length(strs)>1)
    # Guarda a la variable numèrica buida el nº de reads de l'haplotip avaluat 
	  nr[j] <- as.numeric(strs[2])
  }
  # Agrupa els resultats del bucle for en un data frame de 2 columnes i tantes files
  # com haplotips avaluats d'aquella mostra
  IDs <- data.frame(ID=nms,nseqs=nr,stringsAsFactors=FALSE)
  # Guarda el nº de files del data frame (nº haplotips coincidents)
  nall <- nrow(IDs)
  # Guarda el 'Total Number of Reads': sumatori del nº de reads de tots els haplotips
  tnr <- sum(IDs$nseqs)
  # Retorna una llista amb: 
  # 1) La taula de resultats (ID haplotip i nº de reads) per a tots els haplotips
  # 2) Les seqüències dels haplotips
  # 3) El nº total d'haplotips coincidents del fitxer fasta
  # 4) El nº total de reads (sumatori de tots els haplotips)
  return(list(IDs=IDs,seqs=seqs,nall=nall,tnr=tnr))
}

### Funció per ordenar els haplotips en funció del nº de mutacions amb la màster
  # Els arguments són les seqs filtrades per min.pct (en aquest cas 0) i el seu
  # nº de reads
SortByMutations <- function(bseqs,nr)
{ # Indica la seqüència màster, la que presenta major nº de reads dins la mostra
  master <- bseqs[which.max(nr)]
  ## Determinar diferències respecte la màster
  # Agafa les seqs dels haplotips i les transforma a un objecte DNAStringSet
  bseqs <- DNAStringSet(bseqs)
  # Realitza un aliniament global Needleman-Wunsch entre les seqs dels haplotips i 
  # la màster
  psa <- pairwiseAlignment(pattern=bseqs,subject=master)
  # Guarda el nº de mismatches de tots els aliniaments
  nm <- nmismatch(psa)
  # Aplica la funció 'table()' sobre les mutacions per veure quantes seqs presenten
  # el mateix nº de mutacions respecte la màster
  tnm <- table(nm)
  
  ## Ordenar per nombre de mutacions
  # Ordena el nombre de mutacions respecte la màster en ordre ascendent
  o <- order(nm)
  # Ordena les seqs segons el seu nº de mutacions respecte la màster
  bseqs <- bseqs[o]
  # Ordena també el nº de reads segons el nº de mutacions de la seqüència
  nr <- nr[o]
  # Ordena la variable amb el nº de mismatches segons les mutacions (ordre ascendent)
  nm <- nm[o]
  
  ##  Numero d'ordre dins de cada nombre de mutacions
  # length(tnm) calcula el nº d'agrupacions de la taula tnm, és a dir el nº màxim de mutacions que s'han trobat
  # 1:tnm[i] s'aplica sobre els valors de 1 fins al total de mutacions trobades. Per cada nº de mutacions, retorna un 
  # conjunt de nombres que van de l'1 al total de vegades que s'ha donat aquell nº de mutacions 
  # 'unlist()' concatena tots els valors per tots el nº de mutacions
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  
  ## Ordena per freqüència descendent dins de cada nombre de mutacions:
  # as.integer(names(tnm)) retorna els valors de 1 fins al total de mutacions trobades
  for(i in as.integer(names(tnm)))
  { # Indexs de les seqüències que presenten aquell nº de mutacions (i) respecte la màster
    idx <- which(nm==i)
    # Pels valors de freqüència pertanyents a les seqs amb 'i' mutacions, retorna la seva posició
    # assignada en ordenar-les en ordre descendent
    o <- order(nr[idx],decreasing=TRUE)
    # Ordena les seqüències amb 'i' mutacions segons freq en ordre descendent
    bseqs[idx] <- bseqs[idx[o]]
    # Ordena les freqüències de les seqs amb 'i' mutacions en ordre descendent
    nr[idx] <- nr[idx[o]]
  }
  ## Calcula la freqüència relativa dels haplotips
  frq <- round(nr/sum(nr)*100,2)
  
  ## Nom complet per cada haplotip
  # Defineix el nom que rebrà cadascun dels haplotips coincidents en les cadenes
  # nm= mutacions respecte seq màster
  # La funció 'sprintf()' permet concatenar el nº d'ordenació de l'haplotip dins del conjunt
  # amb el mateix nº de mutacions, amb format 0000 (per això "%04d")
  nms <- paste("Hpl",nm,sprintf("%04d",isq),sep="_")
  
  ## Capçalera fasta amb nom, nombre de reads i mutacions
  # Assigna els noms generats dels haplotips a les respectives seqüències
  names(bseqs) <- nms
  # Retorna una llista amb les seqs, el nº de reads i el nº de mutacions
  list(bseqs=as.character(bseqs),nr=nr,nm=nm)
}

### Llegeix les seqüències alineades dels haplotips i retorna el nº de reads
  # i mutacions després d'aliniar amb la màster i tornar a ordenar les seqs
# Els arguments son el fitxer avaluat de la carpeta MACH i el % mínim de reads
# (en aquest script és 0 perquè ja s'han filtrat les seqs)
get.data <- function(flnm,min.pct=0.1) # flnm=file.path(mach.Dir,flnms[i]); min.pct=0
{ # Aplica la funció d'abans per obtenir els haplotips del fitxer i el seu nº de reads
  lst <- read.ampl.seqs(flnm)
  # Indica quins haplotips del fitxer presenten un % de reads major al mínim inicat
  # Nota: sum(lst$IDs$nseqs) és equivalent al valor tnr de la llista (lst$tnr)
  # Nota: en aquest fitxer només es fa servir un cop i min.pct és 0, es podria 
  # eliminar aquest filtre?
  fl <- lst$IDs$nseqs/sum(lst$IDs$nseqs)*100 >= min.pct
  # Actualitza les seqüències dels haplotips amb els que tinguin el nº de reads major
  # a l'indicat
  seqs <- lst$seqs[fl]
  # Els noms de les seqüències correspondran als identificadors dels haplotips,
  # després d'eliminar els que tenen un nº de reads menor
  names(seqs) <- lst$IDs$ID[fl]
  # Actualitza el nº de reads dels haplotips que tinguin el nº de reads major
  # a l'indicat
  nr <- lst$IDs$nseqs[fl]
  # Aplica la funció anterior per obtenir els haplotips ordenats en funció de les
  # mutacions respecte la màster
  lst <- SortByMutations(seqs,nr)
  # Retorna una llista amb les seqs, el nº de reads i el nº de mutacions
  # (mateix resultat que la funció SortByMutations)
  list(seqs=lst$bseqs,nr=lst$nr,nm=lst$nm)
}

#------------------------------------------------------------------------#

### Llegeix el fitxer amb els primers
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                      stringsAsFactors=F)
# Els noms de les files de la taula primers correspondran als pools
rownames(primers) <- primers$Ampl.Nm

## Guarda en una variable el nom dels fitxers de la carpeta MACH que corresponen als
# resultats després de la intersecció dels haplotips forward i reverse 
# (amb l'estructura MACHpl02.EV1.A1.fna). El símbol '+' indica que el caràcter que
# el precedeix ha d'aparèixer com a mínim un cop, i el símbol '$' indica el final
# de la cadena de caràcter a buscar
fl.patt <- "MACHpl02\\..+\\.fna$"

## Llista dels fitxers de la carpeta MACH que coincideixen amb el caràcter
# indicat anteriorment
flnms <- list.files(path=mach.Dir,patt=fl.patt)
# Separa el nom dels fitxers per punts i guarda els resultats en una taula, on 
# les columnes inclouran cada element del nom del fitxer
parts <- t(sapply(flnms,function(str) strsplit(str,split="\\.")[[1]]))
# Els elements de la columna 3 (que abans corresponien al terme HBV) ara seran
# la concatenació de HBV i les coordenades de la regió avaluada (elements de
# les columnes 4 i 5)
parts[,3] <- apply(parts,1,function(x) paste(x[3],x[4],x[5],sep="."))
# Guarda la columna 3 de la taula (concatenació HBV.coordenada1.coordenada2)
anms <- parts[,3]
# Guarda la columna 2 de la taula (ID del pacient)
snms <- parts[,2]

# Guarda el total de fitxers trobats a la carpeta MACH (nº de mostres, 2 per pacient)
n <- length(flnms)
# Genera un data frame amb les columnes indicades a l'argument. La 1a i 2a columna
# inclouran els IDs dels pacients i el nom dels amplicons o regions avaluades, respectivament
rprt <- data.frame(ID=snms,Ampl=anms,Hpl=integer(n),MaxMut=integer(n),
                   Master=numeric(n),Reads=numeric(n),stringsAsFactors=FALSE)

# Bucle sobre el total de mostres (2 per pacient)
for(i in 1:n)
{ # Sobre el fitxer .fna de la mostra avaluada, aplica la funció definida a l'inici 
  # de l'script
  lst <- get.data(file.path(mach.Dir,flnms[i]),0)   
  ## A la taula de resultats rprt guarda les dades:
  # Nº d'haplotips (longitud de la llista del nº de reads)
  rprt$Hpl[i] <- length(lst$nr)
  # Màxim nº de mutacions dels haplotips
  rprt$MaxMut[i] <- max(lst$nm)
  # Calcula el % que representa la seq màster (en nº reads) respecte el total
  # Nota: es pot obtenir de la variable frq de les funcions d'abans
  rprt$Master[i] <- round(max(lst$nr)/sum(lst$nr)*100,2)
  # Total de reads de tots els haplotips
  rprt$Reads[i] <- sum(lst$nr)
}

# Retorna els indexs de la taula rprt derivats d'ordenar les mostres en funció 
# de l'amplicó (regió) avaluat i l'ID del pacient
o <- order(rprt$Ampl,rprt$ID)

## Genera el fitxer .txt que es guardarà a la carpeta de resultats
sink(file.path(resDir,"NtFastasSummary.txt"))
# Inclou la taula de resultats rprt ordenada en funció dels pacients i regions avaluades
cat("\nFasta files summary\n\n")
print(rprt[o,])

# Retorna el resum de les dades de la columna amb el nº de reads per totes les mostres
cat("\nFinal coverage, global summary\n\n")
print(summary(rprt$Reads))

# Retorna el resum de les dades dels reads en funció de l'amplicó o regió avaluat
cat("\nFinal coverage, by amplicon summary\n\n")
print(tapply(rprt$Reads,rprt$Ampl,summary))
sink()

## Genera la paleta de colors pels gràfics
library(RColorBrewer)
pal <- brewer.pal(8,"Dark2")

## Genera el fitxer .pdf amb els gràfics que es guardarà a la carpeta de resultats
pdf(file.path(resDir,"FinalCovBoxplots.pdf"),paper="a4",width=5,height=10)
par(mfrow=c(2,1))

# Calcula el logaritme en base 10 del nº de reads per cada mostra avaluada
y <- log10(rprt$Reads)
# Genera un gràfic boxplot on es representa el log dels reads en funció del pool
# per mostrar la cobertura de cadascun
boxplot(y~rprt$Ampl,border="gray",ylab="Log10 #reads",outline=FALSE,
        ylim=range(y),las=2)
points(jitter(as.integer(factor(rprt$Ampl)),a=0.2),y,pch="+",cex=0.8)	
title(main="Final coverage by amplicon")	
abline(h=log10(10000),col="red",lty=4)

## Retorna els indexs de la taula rprt derivats d'ordenar la columna de reads en
# ordre ascendent
o <- order(rprt$Reads)
# Ordena el nº dels reads en funció dels indexs obtinguts i calcula el logaritme
y <- log10(rprt$Reads[o])
# Fa la coerció a variable factorial de la taula corresponent al nom dels amplicons,
# ordenats segons el nº de reads
g <- factor(rprt$Ampl[o])
# Assigna els colors del gràfic en funció dels pools o amplicons
icol <- pal[as.integer(g)]
# Genera un gràfic de línies verticals indicant el logaritme del nº de reads 
# amb distinció de color segons el pool
plot(1:length(y),y,type="h",col=icol,ylab="log10 #reads",xlab="sorted samples")
legend("topleft",lwd=2,col=pal,legend=levels(g),cex=0.8)
abline(h=log10(10000),col="red",lty=4)

## Genera un altre boxplot per representar el nº de reads segons el pool, aquest
# cop sense calcular el logaritme
boxplot(rprt$Reads~rprt$Ampl,border="gray",ylab="#reads",outline=FALSE,
        ylim=range(rprt$Reads))
points(jitter(as.integer(factor(rprt$Ampl)),a=0.2),rprt$Reads,pch="+",cex=0.8)		
abline(h=10000,col="red",lty=4)

## Torna a calcular el logaritme del nº de reads (sense ordenar) i fa el mateix
# gràfic boxplot del principi, però per alguna raó cada cop que es representa
# surten gràfics diferents(?
y <- log10(rprt$Reads)
boxplot(y~rprt$Ampl,border="gray",ylab="Log10 #reads",outline=FALSE,
        ylim=range(y),las=2)
points(jitter(as.integer(factor(rprt$Ampl)),a=0.2),y,pch="+",cex=0.8)	
title(main="Final coverage by amplicon")	
abline(h=log10(10000),col="red",lty=4)

dev.off()
