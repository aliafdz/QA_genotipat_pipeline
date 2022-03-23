
library(stringr)
library(Biostrings)
#source("./R/seqanalfns.v4.5.R") # D'aquest fitxer només es crida la funció: 

### Funció per llegir les seqüències aliniades dels amplicons
  # Els arguments corresponen al nom del fitxer on es troben les seqs i el valor
  # de mínim nº de reads assignat (per defecte 2)
read.ampl.seqs <- function(flnm,mnr=2) # flnm=flnms[i],mnr=1
{ # Llegeix el fitxer fasta de l'argument i el transforma en un objecte DNAStringSet
  seqs <- as.character(readDNAStringSet(flnm))
  # Guarda els noms de les seqüències del fitxer, és a dir dels haplotips
  # Aquests haplotips són els que han fet intersecció entre ambdues cadenes
  IDstr <- names(seqs)
  # Guarda el total de seqs del fitxer (nº d'haplotips coincidents)
  n <- length(IDstr)
  
  # Variable de caràcter buida amb tantes entrades com haplotips trobi al fitxer
  nms <- character(length=n)
  # Matriu buida amb tantes files com haplotips i 2 columnes, definides amb 'colnames()'
  sts <- matrix(0,nrow=n,ncol=2)
  colnames(sts) <- c("nseqs","pct1") # nº de seqüències i percentatge
  # Bucle sobre els haplotips del fitxer (coincidents)
  for(j in 1:n)
  { # Separa el nom de l'haplotip avaluat segons el símbol "|"
    strs <- strsplit(IDstr[j],split="\\|")[[1]]
    # Guarda a la variable buida d'abans el primer element del nom separat,
    # que correspon a la consecució de "Hpl", el nº de mutacions amb la seq màster
    # i el nº d'ordre dins del grup amb les mateixes mutacions
    nms[j] <- strs[1]
    # Guarda a la matriu buida els dos elements restants del nom de l'haplotip,
    # que corresponen al nº de reads i la freqüència relativa de l'haplotip
    sts[j,] <- as.numeric(strs[2:3])
  }
  # Agrupa els resultats del bucle for en un data frame de 3 columnes i tantes files
  # com haplotips avaluats d'aquella mostra
  IDs <- data.frame(ID=nms,sts,stringsAsFactors=FALSE)
  # Guarda el nº de files del data frame (nº haplotips coincidents)
  nall <- nrow(IDs)
  # Guarda el 'Total Number of Reads': sumatori del nº de reads de tots els haplotips
  tnr <- sum(IDs$nseqs)
  
  ### Filtra segons el nombre mínim de reads (mnr) per haplotip 
  # Indica quines entrades del data frame generat tenen més reads del mínim permès
  flags <- IDs$nseqs >= mnr
  # Retorna una llista amb: 
  # 1) La taula de resultats (ID haplotip, nº de reads i freq) filtrada per mnr
  # 2) Les seqüències dels haplotips filtrats per mnr
  # 3) El nº total d'haplotips coincidents del fitxer fasta
  # 4) El nº total de reads (sumatori de tots els haplotips)
  return(list(IDs=IDs[flags,],seqs=seqs[flags],nall=nall,tnr=tnr))
}

####-------------------------------------------------------------####

### Fitxers a tractar i descripció
# Carrega el fitxer .RData generat en 'RAV_PrimersSplitFromMIDsSplits_WithGaps-pl-v2.30.R'
# (on s'han eliminat els primers de totes les seqüències les quals es classifiquen
# en forward o reverse )
load(file=file.path(repDir,"SplittedReadsFileTable.RData")) # Inclou les taules FlTbl i PoolTbl
# Retorna els indexs de la taula FlTbl derivats d'ordenar les mostres en funció de l'ID del
# pacient, de l'amplicó (regió) avaluat i la cadena forward o reverse 
o <- order(FlTbl$Pat.ID,FlTbl$Ampl.Nm,FlTbl$Str)
# Reordena les entrades de la taula en funció dels indexs anteriors, de manera que s'agrupen
# les entrades corresponents al mateix pacient i en segon ordre a la mateixa regió avaluada
FlTbl <- FlTbl[o,]

### Generació dels noms dels fitxers fasta
# Guarda els indexs de la taula FlTbl segons la cadena asignada als reads
idx.fw <- which(FlTbl$Str=="fw")
idx.rv <- which(FlTbl$Str=="rv")

# Guarda els noms dels fitxers fasta inclosos a la carpeta MACH generats després 
# de fer la intersecció entre cadenes (MACHpl02). 
flnms <- paste("MACHpl02",FlTbl$Pat.ID[idx.fw],
               FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
flnms <- file.path(mach.Dir,flnms)

# Guarda les concatenacions dels ID dels pacients amb la regió amplificada units per guionet
pnms <- paste(FlTbl$Pat.ID[idx.fw],FlTbl$Ampl.Nm[idx.fw],sep=" - ")

## Genera el fitxer .pdf on es guardaran els gràfics de resultats d'aquest codi
 # Es generara un gràfic per cada mostra avaluada
pdf(file=file.path(repDir,"GapsBarPlots.pdf"),paper="a4",width=7,
    height=10)
par(mfrow=c(4,1))
par(mar=c(4.5, 4, 3, 4) + 0.1)

## Bucle sobre totes les mostres (2 per pacient), que coincideix amb el nº de fitxers MACHpl02
for(i in 1:length(flnms))
{ # Si el fitxer de la mostra avaluada no existeix, passa a la següent iteració
  if(!file.exists(flnms[i])) next
  ## Aplica la funció definida (del fitxer, però es pot copiar aquí) per llegir 
  # els fitxers fasta de cada mostra avaluada 
  lst <- read.ampl.seqs(flnms[i],mnr=1) ## Funció de l'arxiu seqanalfns.v4.5.R
  # Guarda l'index de l'haplotip amb màxim nº de reads
  imast <- which.max(lst$IDs$nseqs)
  # Guarda la seqüència de l'haplotip amb màxim nº de reads (màster)
  rsq <- lst$seqs[imast]
  # Guarda la longitud de la seq de l'haplotip màster
  xmx <- nchar(rsq)
  # Guarda totes les seqüències dels haplotips del fitxer avaluat
  seqs <- lst$seqs
  # Guarda la llista amb el nº de reads de tots els haplotips
  nr <- lst$IDs$nseqs
  
  ## Separa la seqüència de l'haplotip màster per nucleòtids 
  rnt <- strsplit(rsq,split="")[[1]]
  # Guarda una variable buida amb tants 0 com nucleòtids hi ha a la seq
  tp <- integer(length(rnt))
  # Indica quins elements de la seq corresponen a insercions (símbol "-")
  ins.pos <- rnt=="-"
  # A la variable buida d'abans, indica amb un 1 les posicions amb insercions
  tp[ins.pos] <- 1  

  ## Ara separa per nucleòtids les seqs de tots els haplotips de la mostra avaluada
  ntmat <- t(sapply(seqs,function(s) strsplit(s,split="")[[1]]))
  # Indica, sobre cada element de les seqs, quins corresponen a delecions ("-")
  del.pos <- apply(ntmat,2,function(vnt) any(vnt=="-"))
  # A la variable d'abans, indica amb un 2 les posicions on no hi ha insercions
  # però sí delecions
  tp[!ins.pos & del.pos] <- 2
  # Suma una unitat a tots els elemnts de la variable tp, de manera que un 1 indica
  # on hi ha nucleòtids, 2 on hi ha insercions i 3 on hi ha delecions
  tp <- tp+1
  
  # Guarda una altra variable buida amb tants 0 com nucleòtids hi ha a la seq màster
  nb <- integer(length(rnt))
  # En cas de que es trobin insercions en alguna de les seqs,
  if( sum(ins.pos) )
  { # A les posicions on s'hagin detectat insercions, en la variable buida nb 
    # es guarda el nº de reads dels haplotips on s'ha donat la inserció (i per tant no hi ha gap)
    nb[ins.pos] <- apply(ntmat[,ins.pos,drop=FALSE],2,function(vnt)
                           sum(nr[vnt!="-"]))
  }
  # En cas de que es trobin delecions en alguna de les seqs,
  if( sum(del.pos) )  
  { # A les posicions on s'hagin detectat delecions i no insercions, en la variable nb 
    # es guarda el total de reads dels haplotips on s'ha donat la deleció (i per tant hi ha gap)
    nb[!ins.pos & del.pos] <- apply(ntmat[,!ins.pos & del.pos,drop=FALSE],2,
                 function(vnt) sum(nr[vnt=="-"]))
  }
  # Genera un gràfic amb la funció 'plot()' per representar els InDels de les seqs dels 
  # haplotips (color vermell insercions i blau delecions), en forma de línies verticals,
  # indicant la posició de la seq en l'eix X i el nº de reads en l'eix Y
  plot(nb,type="h",col=c("black","red","blue")[tp],lwd=2,yaxt="n",
       xlab="MA position",ylab="reads",xlim=c(0,xmx),ylim=c(0,max(nb[1:xmx])))
  # Genera els intervals de l'eix Y de manera que siguin de 20 en 20
  ytk <- axTicks(2)
  # Genera l'eix Y de l'esquerra (nº reads)
  axis(side=2,at=ytk,las=2)
  # Per generar l'eix Y de la dreta (%), es divideixen els valors de l'eix Y 
  # de l'esquerra i es divideix entre el total de reads de la mostra avaluada
  pct <- round(ytk/sum(nr)*100,2)
  axis(side=4,at=ytk,labels=pct,las=2,cex.axis=0.8)
  mtext("Percentage",side=4,line=3,cex=0.6)
  # El títol de cada gràfic correspon a l'ID del pacient i la regió amplificada 
  title(main=pnms[i],line=1)
}

dev.off()
