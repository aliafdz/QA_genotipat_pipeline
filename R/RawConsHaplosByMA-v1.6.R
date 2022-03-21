############################################################
###    RAW CONSENSUS HAPLOTYPES by MULTIPLE ALIGNMENT    ###
############################################################

### Funció per assignar a cada haplotip el seu nom en funció del nº de mutacions que presenta
  # respecte la seqüència màster (més freqüent) i el seu ordre. El nº d'ordenació
  # en tots els casos serà un número de 4 dígits (per això s'afegeixen 0s a l'esquerra)
zeroFillInt2Char <- function(x,ln)
{ # Concatenació de 0 seguits de l'ordre que presenta l'haplotip dins del conjunt
  # d'haplotips amb el mateix nº de mutacions amb la màster
  x <- paste("000000",x,sep="")
  # 'nchar()' recompta el nº de caràcters de la concatenació anterior
  # En aquest cas ln=4. 'substring()' extrau els caràcters de la concatenació
  # des de la posició indicada (inclosa) fins al final. Així assegura l'addició de 0s
  # fins que el nombre total tingui 4 dígits
  substring(x,nchar(x)-ln+1)
}


### Funció per separar els noms dels haplotips amb el recompte de reads i la freqüència
split.fasta.names <- function(seqs,mnr=2)
{ # Guarda els noms assignats a les seqüències de l'objecte DNAStringSet de l'argument
  IDstr <- names(seqs)
  # Guarda el total de noms assignats del pas anterior, és a dir el total d'haplotips
  n <- length(IDstr)
  # Guarda un vector buit amb tantes entrades com haplotips
  nms <- character(length=n)
  # Guarda una matriu buida amb n files (1 per haplotip) i dues columnes, definides amb 
  # 'colnames()', on es guardaran en nº de reads de l'haplotip i la seva freqüència
  sts <- matrix(0,nrow=n,ncol=2)
  colnames(sts) <- c("nseqs","pct1")
  # Iteració sobre el total d'haplotips del fitxer avaluat
  for(j in 1:n)
  { # Separa el nom de l'identificador de l'haplotip, que inclou 3 elements:
    # 1) Nº de mutacions i ordre de l'haplotip dins del grup amb aquestes mutacions
    # 2) Nº de reads de l'haplotip
    # 3) Freq relativa dins de la mostra
    strs <- strsplit(IDstr[j],split="\\|")[[1]]
    # Guarda a la taula nms el primer element de l'identificador de l'haplotip corresponent
    nms[j] <- strs[1]
    # Guarda a la 1a i 2a columna de la taula sts el 2n i 3r elements de l'identificador
    sts[j,] <- as.numeric(strs[2:3])
  }
  # Guarda en un data frame els 3 elements de l'identificador de tots els haplotips separats per
  # columnes
  # NOTA: es podria simplificar per tal de guardar les variables directament a partir del bucle?
  IDs <- data.frame(ID=nms,sts,stringsAsFactors=FALSE)
  # Guarda el nº de files del data frame (total d'haplotips)
  nall <- nrow(IDs)
  # Realitza el sumatori de la columna amb el nº de reads per haplotip
  tnr <- sum(IDs$nseqs)
  # Guarda el nº de reads majors al valor definit de minim nombre de reads
  flags <- IDs$nseqs >= mnr
  # Retorna una llista amb un subset dels haplotips que presentin un nº de reads major al mínim definit,
  # les seqüències dels haplotips que compleixen la condició, el nº de total d'haplotips i el sumatori de reads
  return(list(IDs=IDs[flags,],seqs=seqs[flags],nall=nall,tnr=tnr))
}


### Funció per llegir els .fna amb els amplicons ja classificats i sense primers 
read.ampl.seqs <- function(flnm,mnr=2)
{ # Guarda les seqüències del fitxer fasta de l'argument obtingudes amb la funció 'readDNAStringSet()'
  seqs <- as.character(readDNAStringSet(flnm))
  # Aplica la funció definida anteriorment per obtenir tots els noms dels haplotips amb un nº de reads
  # major al mínim definit i les seves dades
  return(split.fasta.names(seqs,mnr)) # NOTA: es podrien fusionar les 2 funcions?
}


###  Escriu unes seqüències a un fitxer fasta
############################################### NOTA: aquesta i la següent es poden eliminar
write.fasta <- function(seqs,flnm)            # i fer servir la funció de R
{ writeXStringSet(seqs,flnm) }

###  Llegeix les seqüències d'un fitxer fasta
############################################### 
read.fasta <- function(flnm)
{ readDNAStringSet(flnm) }


###  Funció per a l'aliniament múltiple amb muscle de les seqs
###  Les seqüècies s'entren i es tornen com a DNAStringSet
# Guarda la ruta de l'executable muscle (aquest es posarà al fitxer global)
muscle <- "C:\\Muscle\\muscle3.8.31_i86win32.exe"
# Guarda el nom del fitxer d'opcions de muscle que es guardarà a l'entorn global del projecte
muscle.cl.opts <- c("-log muscle.log")

doMuscle <- function(seqs)
{ # Genera un fitxer al directori tmp per introduir les seqüències a alinear
  tmp.file <- file.path(tmp.Dir,"muscleInFile.fna")
  # Genera el fitxer al directori tmp que resultarà de l'alineament amb muscle
  res.file <- file.path(tmp.Dir,"muscleOutFile.fna")
  # Si ja existeix l'arxiu de resultats, s'elimina per generar-lo de nou
  if(file.exists(res.file)) file.remove(res.file)
  # Aplica la funció d'abans per escriure el fitxer fasta d'entrada amb les seqüències de l'argument
  write.fasta( seqs, tmp.file )
  # Genera els fitxers in i out per a executar muscle
  in.file <- paste("-in ",tmp.file,sep="")
  out.file <- paste("-out ",res.file,sep="")
  # Genera la comanda per executar muscle amb el fitxer .exe, els arxius in i out i el fitxer d'opcions
  command <- paste(muscle,in.file,out.file,
                paste(muscle.cl.opts,collapse=" "),sep=" ")
  # Empra la funció 'system()' per executar la comanda indicada abans i fer l'aliniament múltiple
  system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
          ignore.stdout=FALSE,invisible=TRUE)
  # Retorna el fitxer de resultats com a objecte DNAStringSet gràcies a la funció de llegir arxius fasta
  if( file.exists(res.file) )
    return( read.fasta(res.file) )
  return(NULL)
}


###  Taula de freqüències de nucleòtid per posició
#####################################################
SeqsPosTable <- function(seqs,nr)
{ col.freqs<- function(s,nr) 
  { c( sum(nr[s=="A"]),sum(nr[s=="C"]),sum(nr[s=="G"]),
       sum(nr[s=="T"]),sum(nr[s=="-"]) )
  }
  nuc.mat <- t(sapply(seqs,function(s) strsplit(s,split="")[[1]]))
  frpos <- t(apply(nuc.mat,2,col.freqs,nr))
  colnames(frpos) <- c("A","C","G","T","-")
  frpos
}

       
### Alinear les distribucions d'haplotips de dues poblacions
# Els arguments de la funció corresponen al nº de seqüències dels haplotips FW o RV (nA o nB)
# i les seqüències en si (seqsA i seqsB)
PopsAlgnHist <- function(nA,seqsA,nB,seqsB)
{ # Relaciona el nom dels haplotips amb el nº total de seqüències que inclouen
  names(nA) <- seqsA
  names(nB) <- seqsB
  # Combina les dades dels dos conjunts de seqüències, eliminant aquelles que estiguin duplicades en algun conjunt
  nms <- union(seqsA,seqsB) # tots: A + B
  # Guarda la longitud (nº de seqs) de la nova combinació -> No correspon a la suma de seqsA i seqsB perquè elimina els duplicats
  nb <- length(nms)
  # Guarda un vector buit de longitud igual al nº de seqüències totals, per guardar les coincidències amb FW
  pA <- integer(nb)
  # Busca quines seqüències del conjunt global coincideixen amb els haplotips FW
  idx <- which(nms %in% seqsA)
  # Guarda en el vector d'abans el nº de seqüències dels haplotips que han coincidit en la seva respectiva posició
  pA[idx] <- nA[nms[idx]]
  # Guarda un vector buit de longitud igual al nº de seqüències totals, per guardar les coincidències amb RV
  pB <- integer(nb)
  # Busca quines seqüències del conjunt global coincideixen amb els haplotips RV
  idx <- which(nms %in% seqsB) 
  # Guarda en el vector d'abans el nº de seqüències dels haplotips que han coincidit en la seva respectiva posició
  pB[idx] <- nB[nms[idx]]
  # Genera una llista amb el nº de seqs dels haplotips FW i RV, i també el conjunt global amb totes les seqs
  list(pA=pA,pB=pB,Hpl=nms)
}


### Confrontació d'haplotips FW i RV
# Els arguments de la funció corresponen al nº de seqüències dels haplotips FW o RV (nA o nB)
# i les seqüències en si (seqsA i seqsB)
Intersect.FWRV <- function(nA,seqsA,nB,seqsB)  
{ ### Aliniament de distribucions
  # Aplica la funció d'abans per obtenir el conjunt de totes les seqs d'haplotips FW i RV i el seu nº de reads 
  lst <- PopsAlgnHist(nA,seqsA,nB,seqsB)
  
  ### Haplotips comuns i solapament global
  # Indica en quins casos coincideixen les seqs dels haplotips FW i RV, ja que seran els indexs on el nº de reads
  # per A i per B sigui major a 0
  fl <- lst$pA>0 & lst$pB>0
  # Realitza el sumatori de seqs dels haplotips coincidents en FW i RV i el divideix entre el total
  ov.a <- sum(lst$pA[fl]+lst$pB[fl])/sum(lst$pA+lst$pB)
  
  ### Freqüències i solapament per intersecció
  # Calcula la freq relativa (nº de seqs d'un haplotip entre el total) per aquells haplotips coincidents en ambdues cadenes
  # Es considera que el sumatori de reads en ambdues cadenes correspon a la cobertura de l'amplicó
  pFW <- (lst$pA/sum(lst$pA))[fl]
  pRV <- (lst$pB/sum(lst$pB))[fl]
  # Indica el valor mínim de freq al comparar cada parella d'haplotips coincidents -> intersecció
  p <- pmin(pFW,pRV) 
  # Calcula el sumatori de tots els valors mínims de freq calculats entre les parelles coincidents
  ov.i <- sum(p)
  # Ara renormalitza les dades dividint cada valor mínim de freq entre el sumatori de totes 
  p <- p/sum(p)
  # Retorna una llista que inclou, en ordre:
  # 1) Freq relatives dels haplotips normalitzades
  # 2) Seqüències dels haplotips coincidents en FW i RV
  # 3) Valor sumatori de tots els valors mínims de freq calculats entre les parelles coincidents
  # 4) Valor sumatori del nº de seqs dels haplotips coincidents en FW i RV entre el total
  # 5) Vector amb el nº de seqs dels haplotips FW presents en el conjunt global
  # 6) Vector amb el nº de seqs dels haplotips RV presents en el conjunt global
  list(p=p,seqs=lst$Hpl[fl],ov.i=ov.i,ov.a=ov.a,pA=lst$pA,pB=lst$pB)
}


### Funció per representar les distribucions aliniades dels haplotips
PlotHplHistos <- function(tt,pA,pB,p)
{ # Divideix el nº de seqs de cada haplotip de la cadena FW entre el total
  pA <- pA/sum(pA)
  # Divideix el nº de seqs de cada haplotip de la cadena RV entre el total
  pB <- pB/sum(pB)
  # Defineix el límit de l'eix Y en funció de les freqüències calculades
  ymx <- max(c(pA,pB))
  # Gràfic de barres per representar els haplotips de la cadena FW i la seva freq relativa
  barplot(pA,ylim=c(0,ymx)); abline(h=0)
  title(main="FW strand haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  # Gràfic de barres per representar els haplotips de la cadena RV i la seva freq relativa
  barplot(pB,ylim=c(0,ymx)); abline(h=0)
  title(main="RV strand haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
  # Gràfic de barres per representar els haplotips coincidents en ambdues cadenes i la seva freq
  barplot(p); abline(h=0)
  title(main="Intersected haplotypes barplot",line=0.5,cex.main=1)
  title(sub=tt,line=0.5,cex.main=1)
}



### Dóna nom als haplotips en funció del nombre de mutacions respecte la màster
### i de la seva frequencia poblacional, i salva les seqüències en un fitxer fasta.
# Els arguments corresponen al nom del fitxer on es desaran els resultats, les seqs
# dels haplotips coincidents en ambdues cadenes, el sumatori de reads per haplotip coincident
# i el nº màxim de diferències permeses
flnm=out.flnms[i]
bseqs=lst$seqs
nr=rds
SaveHaplotypes <- function(flnm,bseqs,nr,max.difs=250)  #,seq0)
{
  ## Si només hi ha un haplotip coincident:
  if(length(bseqs)==1)
  { # Separa la seqüència per nucleòtids
    bnts <- strsplit(bseqs[1],split="")[[1]]
    # Guarda únicament els nucleòtids de la seqüència, elimina els gaps (-) 
  	bnts <- bnts[bnts!="-"]
  	# Torna a generar la seqüència de DNA tornant a unir els nucleòtids
  	bseqs <- paste(bnts,collapse="")
  	# Assigna a la seq el nom corresponent: haplotip nº1 amb 0 mutacions, indicant el nº de reads, 
  	# i la freq serà del 100% (perquè només queda 1 haplotip)
  	names(bseqs) <- paste("Hpl.0.0001",nr,100,sep="|")
  	# Guarda la seq en un fitxer fasta en la ruta indicada en l'argument de la funció
    write.fasta(DNAStringSet(bseqs),flnm)
  # Retorna una llista amb la seqüència i el nº de reads de l'haplotip
	return( list(bseqs=bseqs,nr=nr,nm=0) )
  }

  ## Determinar diferències respecte la seqüència (haplotip) màster
  # bseqs[which.max(nr)] permet assignar la màster com l'haplotip amb major nº de reads
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[which.max(nr)])
  # Guarda el nº de diferències obtingudes dels aliniaments
  nm <- nmismatch(psa)

  ## Elimina les seqs amb massa diferències respecte la màster
  # Subset de seqüències d'haplotip que tenen un nombre de diferències menor al màxim permès
  bseqs <- bseqs[nm<=max.difs]
  # També actualitza el nº de reads només dels haplotips amb menys diferències del màxim permès
  nr <- nr[nm<=max.difs]
  # Elimina les dades dels aliniaments amb major nº de diferències del permès
  nm <- nm[nm<=max.difs]
  
  # Si després d'eliminar els haplotips amb múltiples diferències només en queda un, realitza el procés
  # del principi per guardar aquella seqüència en el fitxer fasta
  if(length(bseqs)==1)
  { bnts <- strsplit(bseqs[1],split="")[[1]]
	bnts <- bnts[bnts!="-"]
	bseqs <- paste(bnts,collapse="")
	names(bseqs) <- paste("Hpl.0.0001",nr,100,sep="|")
    write.fasta(DNAStringSet(bseqs),flnm)
	return( list(bseqs=bseqs,nr=nr,nm=0) )
  }

  ## Ordenar per nombre de diferències
  # Ordena el nombre de mutacions respecte la màster en ordre ascendent
  o <- order(nm)
  # Ordena les seqs segons el seu nº de mutacions respecte la màster
  bseqs <- bseqs[o]
  # Ordena també les freqüències segons el nº de mutacions de la seqüència
  nr <- nr[o]
  # Ordena la variable amb el nº de mismatches segons les mutacions (ordre ascendent)
  nm <- nm[o]

  ## Numero d'ordre dins de cada nombre de mutacions:
  # Guarda les vegades que apareix cada nº de mismatches
  tnm <- table(nm)
  # length(tnm) calcula el nº d'agrupacions de la taula tnm, és a dir el nº màxim de mutacions que s'han trobat
  # 1:tnm[i] s'aplica sobre els valors de 1 fins al total de mutacions trobades. Per cada nº de mutacions, retorna un 
  # conjunt de nombres que van de l'1 al total de vegades que s'ha donat aquell nº de mutacions 
  # 'unlist()' concatena tots els valors per tots el nº de mutacions
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  
  ## Ordena per freqüència descendent dins de cada nombre de mutacions:
  # as.integer(names(tnm) retorna els valors de 1 fins al total de mutacions trobades
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

  ## Calcula la freqüència relativa dels haplotips coincidents en FW i RV
  frq <- round(nr/sum(nr)*100,2)

  ## Nom complet per cada haplotip
  # Defineix el nom que rebrà cadascun dels haplotips coincidents en les cadenes
  # nm= mutacions respecte seq màster
  # 'zeroFillInt2Char' definida al principi, retorna el nº d'ordenació de l'haplotip dins del conjunt
  # amb el mateix nº de mutacions
  nms <- paste("Hpl",nm,zeroFillInt2Char(isq,4),sep=".")

  ## Genera la capçalera fasta amb el nom de l'haplotip, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ###  Afegir-hr la RefSeq -> no s'executa
  #bseqs <- c(seq0,bseqs)

  ### Eliminar columnes de tot gaps
  # Per cadascuna de les seqs dels haplotips coincidents, aplica la funció 'strsplit()' per
  # separar-les per nucleòtids
  nuc.mat <- t(sapply(bseqs,function(s) strsplit(s,split="")[[1]]))
  # Guarda els gaps de totes les seqüències
  fl <- apply(nuc.mat,2,function(st) all(st=="-"))
  # Si es troba algun gap, 
  if(sum(fl))
  { # Guarda només els nucleòtids sense gaps
    nuc.mat <- nuc.mat[,!fl]
    # Torna a ensamblar els nucleòtids per generar les seqs de DNA
    bseqs <- apply(nuc.mat,1,paste,collapse="")
  }

  ## Genera el fitxer fasta amb els haplotips que han coincidit i sense gaps
  write.fasta(DNAStringSet(bseqs),flnm)
  # Retorna una llista amb les seqs dels haplotips coincidents, el nº de reads i el nº de mutacions
  list(bseqs=bseqs,nr=nr,nm=nm)
}


###  Mesures de diversitat d'un amplicó
########################################
AmplStats <- function(seqs,nr)
{
  smres <- data.frame( Hapl=integer(1),Eta=integer(1),
                       S=integer(1),Mf=numeric(1),Sn=numeric(1),
                       Pi=numeric(1) )
  smres$Hapl <- length(seqs)
  ###  Compute the nucleotide frequencies and the mutations table
  pos.tbl <- SeqsTbl.w(seqs,nr)
  mut.tbl <- MutsTbl(pos.tbl)
  ### Total number of mutations (Eta)
  smres$Eta <- TotalMutations(mut.tbl)
  ### Number of polymorphic sites (segregating sites, S)
  smres$S <- PolymorphicSites(mut.tbl)
  ### Mutations frequency
  smres$Mf <- TotalMutations.w(mut.tbl) / 
                (sum(nr)*nchar(seqs[1]))
  ### Normalized quasispecies Shannon entropy
  smres$Sn <- QSEntropy(nr)/log(smres$Hapl)
  ###  Pairwise differences
  d.ij <- PairwiseDiffs.w(seqs,nr)
  ###  Nucleotide diversity (Pi)
  smres$Pi <- MeanNuclDiffs.w(d.ij,nr)/nchar(seqs[1])
  smres
}

#-------------------------------------------------------------------#

library(Biostrings)
library(ape)

# En aquest fitxer trobem tot tipus de funcions definides per l'anàlisi d'aquest script
source("./R/seqanalfns.v4.5.R")

### Llegeix l'estructura de descripció de mostres
samples <- read.table(file.path(desc.Dir,"samples.csv"), sep="\t",
                      header=T,stringsAsFactors=F)
### Llegeix el fitxer amb els primers
primers <- read.table(file.path(desc.Dir,"primers.csv"), sep="\t",
                      header=T,stringsAsFactors=F)

# Aquest codi no s'executa, es podria borrar?
# RefSeqs <- read.table(file.path(desc.Dir,"RefSeqs.csv"), sep="\t",
                      # header=T,stringsAsFactor=F)

# Carrega el fitxer RData generat en el pas anterior del pipeline, on s'han eliminat els primers
# de totes les seqüències les quals es classifiquen en forward o reverse 
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
# Retorna els indexs de la taula FlTbl derivats d'ordenar les mostres en funció de l'ID del
# pacient, de l'amplicó (regió) avaluat i la cadena forward o reverse 
o <- order(FlTbl$Pat.ID,FlTbl$Ampl.Nm,FlTbl$Str)
# Reordena les entrades de la taula en funció dels indexs anteriors, de manera que s'agrupen
# les entrades corresponents al mateix pacient i en segon ordre a la mateixa regió avaluada
FlTbl <- FlTbl[o,]

### Guarda els noms dels fitxers fasta inclosos al directori trim i els separa en funció
# de la cadena asignada als reads
in.files <- file.path(trimDir,FlTbl$File.Name)
idx.fw <- which(FlTbl$Str=="fw")
idx.rv <- which(FlTbl$Str=="rv")

# Guarda els noms dels fitxers resultants d'aquest script, que correspondran a la concatenació
# del terme MACHpl02, l'ID del pacient i el nom de la regió avaluada, indicats de manera que només
# s'obtingui un fitxer .fna per a cada regió del mateix pacient (per això el condicional fw)
out.flnms <- paste("MACHpl02",FlTbl$Pat.ID[idx.fw],
               FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
# Guarda la ruta on es desaran els fitxers .fna, a la carpeta MACH
out.flnms <- file.path(mach.Dir,out.flnms)

# Guarda els noms d'uns altres fitxers resultants, que correspondran a la concatenació
# del terme MAfwrv, l'ID del pacient i el nom de la regió avaluada, indicats de manera que només
# s'obtingui un fitxer .fna per a cada regió del mateix pacient (per això el condicional fw)
ma.flnms <- paste("MAfwrv",FlTbl$Pat.ID[idx.fw],
               FlTbl$Ampl.Nm[idx.fw],"fna",sep=".")
# Guarda la ruta on es desaran els fitxers .fna, a la carpeta MACH
ma.flnms <- file.path(mach.Dir,ma.flnms)

# Guarda la meitat de la longitud de fitxers presents a la carpeta trim del fitxer previ 
n <- length(in.files)/2

# Genera un data frame inicialment buit, amb n files (2 entrades per pacient) i 5 columnes
# que inclouran els resultats obtinguts en pasos posteriors per a les cadenes fw
rdf.fw <- data.frame(fw.all=integer(n),fw.lowf=integer(n),
                     fw.in=integer(n),fw.unq=integer(n),fw.com=integer(n))
# Genera un altre data frame també buit, amb n files i 5 columnes per incloure resultats
# de les cadenes rv
rdf.rv <- data.frame(rv.all=integer(n),rv.lowf=integer(n),
                     rv.in=integer(n),rv.unq=integer(n),rv.com=integer(n))
# Genera un altre data frame buit, amb n files i 6 columnes
rdf.gbl <- data.frame(all=integer(n),lowf=integer(n),unq=integer(n),
                      ovrlp=numeric(n),common=numeric(n),Fn.rd=integer(n))

# Genera el primer fitxer pdf on s'inclouran diverses representacions de resultats generats en el bucle for
# Per cada regió amplificada de cada pacient, es representen els haplotips de cada cadena i els que han 
# coincidit entre les dues, amb les corresponents freqüències
pdf(file.path(repDir,"MA.Intersects.plots.pdf"),paper="a4",
    width=7,height=11)
par(mfrow=c(3,1))

# Bucle for per iterar sobre totes les mostres (2 per pacient)
for(i in 1:n)
{ # Si no existeix el fitxer de la mostra avaluada a la carpeta trim, torna a començar la iteració següent
  if(!file.exists(in.files[idx.fw[i]])) next
  if(!file.exists(in.files[idx.rv[i]])) next

  ## Filtrat implícit per mnr = min.rd (mínim nombre de reads, definit al fitxer de paràmetres)
  # Aplica la funció definida al principi per obtenir els haplotips del fitxer .fna avaluat,
  # en concret pels haplotips de la cadena forward
  lst1 <- read.ampl.seqs(in.files[idx.fw[i]],mnr=min.rd) 
  # Guarda el vector que inclou el nº de seqüències dels haplotips
  nr1 <- lst1$IDs$nseqs
  # Guarda les seqüències dels haplotips amb més reads del mínim permès
  seqs <- lst1$seqs
  ## Guarda les seqüències de longitud major al mínim permès (definit al fitxer de paràmetres)
  fl <- nchar(seqs) >= min.seq.len
  # Filtra les seqüències per eliminar les que presenten longitud menor al mínim permès
  nr1 <- nr1[fl]
  seqs <- seqs[fl]
  ## Filtrar per mínima abundància
  # a.cut (definit al fitxer principal) correspon al min d'abundancia per entrar 
  # en l'aliniament múltiple (%)
  # Guarda les seqüències que presenten una abundància major al mínim permès
  fl1 <- nr1/sum(nr1)*100 >= a.cut
  ## Guarda a la taula de resultats per a cadenes fw:
  # El total de reads de la mostra avaluada
  rdf.fw$fw.all[i]  <- sum(nr1)
  # El nº de reads amb baixa freqüència (menor al mínim permès)
  rdf.fw$fw.lowf[i] <- sum(nr1[!fl1])
  # Filtra les seqüències per eliminar les que presenten freq menor a la mínima permesa
  seqs <- seqs[fl1]
  # Substitueix amb la funció 'sub()' el terme Hpl per HplFw en els noms de les seqüències
  # de la mostra avaluada
  names(seqs) <- sub("Hpl","HplFw",names(seqs))
  # Renom de la variable
  aseqs <- seqs
  # Guarda la longitud del primer haplotip de la mostra
  rawln <- nchar(seqs[1])

  ## Aplica el mateix procés per a les cadenes classificades reverse de la mostra avaluada
  ## Filtrat implicit per mnr = min.rd
  # Aplica la funció definida al principi per obtenir els haplotips del fitxer .fna avaluat,
  # en aquest cas dels haplotips de la cadena reverse de la mateixa mostra avaluada
  lst2 <- read.ampl.seqs(in.files[idx.rv[i]],mnr=min.rd) 
  nr2 <- lst2$IDs$nseqs
  seqs <- lst2$seqs
  ## Eliminar les seqüències més curtes del mínim permès
  fl <- nchar(seqs) >= min.seq.len
  nr2 <- nr2[fl]
  seqs <- seqs[fl]
  ## Filtrar per mínima abundància
  fl2 <- nr2/sum(nr2)*100 >= a.cut
  ## Afegeix els resultats a la taula per cadenes reverse
  rdf.rv$rv.all[i]  <- sum(nr2)
  rdf.rv$rv.lowf[i] <- sum(nr2[!fl2])
  seqs <- seqs[fl2]
  names(seqs) <- sub("Hpl","HplRv",names(seqs))
  
  ## Afegeix a la variable d'abans (amb les cadenes fw) les cadenes rv filtrades
  aseqs <- c(aseqs,seqs)
  # Guarda l'identificador del primer emprat per a l'amplificació de la mostra avaluada
  ipr <- FlTbl$Pr.ID[idx.fw[i]]

  
  ###  Aliniament múltiple per muscle dels haplotips fw, rv i seqüències RefSeq
  # Guarda el resultat de l'aliniament múltiple realitzat amb la funció definida al principi
  seqs <- doMuscle(DNAStringSet(aseqs))
  # Copia el fitxer resultant de l'aliniament múltiple al directori MACH per guardar l'aliniament 
  # dels haplotips de la mostra avaluada
  file.copy(file.path(tmp.Dir,"muscleOutFile.fna"),ma.flnms[i],
            overwrite=TRUE)

  # Aplica la funció per obtenir els noms dels haplotips i les seves dades
  lst <- split.fasta.names(as.character(seqs))
  # Guarda el vector que inclou el nº de seqüències dels haplotips
  nr <- lst$IDs$nseqs
  # Guarda les seqüències dels haplotips amb més reads del mínim permès (per defecte a la funció mnr=2)
  seqs <- lst$seqs

  ## Separar seq FW i RV ja aliniades
  # Després de l'aliniament múltiple de totes les seqüències de la mostra, separa en 2 variables les que 
  # corresponen a cadena forward o reverse
  ifw <- grep("^HplFw",names(seqs))
  irv <- grep("^HplRv",names(seqs))
  
  # Aplica una altra funció definida al principi per calcular la intersecció entre els haplotips
  lst <- Intersect.FWRV(nr[ifw],seqs[ifw],nr[irv],seqs[irv])
  # Indica en quins casos coincideixen les seqs dels haplotips FW i RV, ja que seran els indexs on el nº de reads
  # per A i per B sigui major a 0
  fl <- lst$pA>0 & lst$pB>0  # Nota: Aquesta variable es podria obtenir directament de la funció

#### Aquest codi no s'executa
  # cat("\n  Haplotypes distribution overlap: ",round(lst$ov.i*100,2),"%",
      # sep="")
  # cat("\n       Overlap as reads in common: ",round(lst$ov.a*100,2),"%\n",
      # sep="")
  # cat("\n       Input reads: ",sum(lst$pA),", ",sum(lst$pB),sep="")
  # cat("\n   Reads in common: ",sum(lst$pA[fl]),", ",sum(lst$pB[fl]),sep="")
  # cat("\n        Reads lost: ",sum(lst$pA[!fl]),", ",
                               # sum(lst$pB[!fl]),sep="")
  # cat("\n\n")

  ## Guarda a les taules de resultats de cadenes FW i RV:
  # El nº total de seqs dels haplotips coincidents en el conjunt global d'haplotips FW+RV
  rdf.fw$fw.in[i]  <- sum(lst$pA)
  rdf.rv$rv.in[i]  <- sum(lst$pB)
  # El nº de seqs dels haplotips coincidents entre ambdues cadenes: cobertura de l'amplicó
  rdf.fw$fw.com[i] <- sum(lst$pA[fl])
  rdf.rv$rv.com[i] <- sum(lst$pB[fl])
  # El nº de seqs dels haplotips no coincidents entre ambdues cadenes
  rdf.fw$fw.unq[i] <- sum(lst$pA[!fl])
  rdf.rv$rv.unq[i] <- sum(lst$pB[!fl])
  
  ## Guarda a la taula de resultats globals (que en aquest punt encara està buida):
  # El total de reads FW+RV de la mostra avaluada
  rdf.gbl$all[i]    <- rdf.fw$fw.all[i]+rdf.rv$rv.all[i]
  # Sumatori de tots els valors mínims de freq calculats entre les parelles coincidents: intersecció entre reads
  rdf.gbl$ovrlp[i]  <- round(lst$ov.i*100,2)
  # Sumatori del nº de seqs dels haplotips coincidents en FW i RV entre el total: superposició haplotips FW i RV
  rdf.gbl$common[i] <- round(lst$ov.a*100,2)
  # Sumatori del nº de seqs totals dels haplotips coincidents d'ambdues cadenes
  rdf.gbl$Fn.rd[i]  <- rdf.fw$fw.com[i]+rdf.rv$rv.com[i]
  # Sumatori d reads dels haplotips amb baixa freqüència d'ambues cadenes
  rdf.gbl$lowf[i]   <- rdf.fw$fw.lowf[i] + rdf.rv$rv.lowf[i]
  # Sumatori dels reads únics d'una cadena de DNA (no coincidents)
  rdf.gbl$unq[i]    <- rdf.fw$fw.unq[i] + rdf.rv$rv.unq[i]

  ### Si no es detecta superposició entre els haplotips d'ambdues cadenes, salta a la seqüent iteració
  if(sum(fl)==0) 
  { cat("\n--------------------------------------------------\n")
    next
  }

  ### Gràfic de barres dels haplotips alineats
  # Concatenació de l'ID del pacient i la regió avaluada en la iteració
  tt <- paste(FlTbl$Pat.ID[idx.fw[i]]," - ",FlTbl$Ampl.Nm[idx.fw[i]],sep="")
  # Suma el nº de seqs dels haplotips FW i RV independentment de si coincideixen o no
  p <- lst$pA+lst$pB
  # Substitueix el nº de seqs d'aquells haplotips no coincidents per 0: Només queda el sumatori de reads
  # dels haplotips que han coincidit en ambdues cadenes
  p[ lst$pA==0 | lst$pB==0 ] <- 0
  # Calcula la freq relativa del nº de seqs o reads de cada haplotip entre el total dels que han coincidit
  p <- p/sum(p)
  # Aplica la funció local per representar els haplotips de cada cadena i els aliniats amb les seves freq
  PlotHplHistos(tt,lst$pA,lst$pB,p)

  ### Suma de reads dels haplotips FW i RV coincidents,
  rds <- lst$pA[fl]+lst$pB[fl]

  ### Aplica la funció local del principi per guardar els haplotips coincidents en ambdues cadenes
  # amb les seves freqüències
  #Save common haplos with estimated frequencies
  SaveHaplotypes(out.flnms[i],lst$seqs,rds)  #,seq0)
}

dev.off()


sink(file=file.path(repDir,"MA.Intersects-SummRprt.txt"))

cat("\n   FW + RV HAPLOTYPES INTERSECTIONS")
cat("\n======================================\n")
cat("\nCutting FW and RV at ",a.cut,"% followed by haplotypes intersection.\n",
    sep="")
cat("\nFrequencies as sum of reads FW+RV.\n\n")

cat("\nSUMMARY RESULTS BY READS NUMBER")
cat("\n===============================\n\n")
fl.fw <- FlTbl$Str=="fw"
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.fw,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.rv,stringsAsFactors=FALSE)
print(frdf)
cat("\n")
frdf <- data.frame(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw],
                  rdf.gbl,stringsAsFactors=FALSE)
print(frdf)

cat("\n\nTotal counts:\n\n")
tots.fw <- apply(rdf.fw,2,sum)
print(tots.fw)
cat("\n")
tots.rv <- apply(rdf.rv,2,sum)
print(tots.rv)
cat("\n")
tots.gbl <- apply(rdf.gbl[,c(1:3,6)],2,sum)
print(tots.gbl)
cat("\n")

cat("\nPercentage:\n\n")
ptts <- round(tots.fw/tots.fw[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.rv/tots.rv[3]*100,2)
print(ptts)
cat("\n")
ptts <- round(tots.gbl/tots.gbl[1]*100,2)
print(ptts)
cat("\n")

sink()

save(frdf,file=file.path(repDir,"IntersectedReads.RData"))

pdf.flnm2 <- "IntersectBarplots.pdf"
pdf(file.path(repDir,pdf.flnm2),paper="a4r",width=10,height=6)
par(mar=c(7,4,4,2)+.1)

pal <- cls[1:2]
bnms <- paste(Pat.ID=FlTbl$Pat.ID[fl.fw],Ampl.Nm=FlTbl$Ampl.Nm[fl.fw])
vals <- cbind(rdf.fw$fw.com,rdf.rv$rv.com)
ymx <- max(rowSums(vals))*1.2
barplot(t(vals),col=pal,ylim=c(0,ymx),names.arg=bnms,
        las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("fw","rv"),cex=0.8)
title(main="Intersected reads")

vals <- cbind(rowSums(rdf.fw[,c("fw.lowf","fw.unq")]),
              rowSums(rdf.rv[,c("rv.lowf","rv.unq")]))
ymx <- max(rowSums(vals))*1.2
barplot(t(vals),col=pal,ylim=c(0,ymx),names.arg=bnms,
        las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("fw","rv"),cex=0.8)
title(main="Filtered out reads")
			  
pal <- cls[c(1,3,2,4)]
mbp <- data.matrix(frdf[,c(8,5,4)])
mbp.nms <- paste(frdf$Pat.ID,frdf$Ampl.Nm,sep=".")
omar <- par(mar=c(6,3.5,3,2))
barplot(t(mbp),col=pal,ylim=c(0,max(rowSums(mbp))*1.2),names.arg=mbp.nms,
    las=2,cex.names=0.6,cex.axis=0.8)
legend("top",horiz=TRUE,fill=pal,legend=c("Common","Unique","Low freq."),
        cex=0.8)
par(mar=omar)
barplot(frdf$common,col="lavender",ylim=c(0,100),names.arg=mbp.nms,
        las=2,cex.names=0.6,cex.axis=0.8)
title(main="FW + RV intersection yield",line=1)

###  Carregar PoolTbl i Fltbl 
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))

all.nms <- paste(FlTbl$Pat.ID,FlTbl$Ampl.Nm,sep=".")
trds <- tapply(FlTbl$Reads,all.nms,sum)
all.nms <- names(trds)
frds <- frdf$Fn.rd
names(frds) <- paste(frdf$Pat.ID,frdf$Ampl.Nm,sep=".")
ytbl <- matrix(0,nrow=length(all.nms),ncol=2)
rownames(ytbl) <- all.nms
colnames(ytbl) <- c("All","Passed")
ytbl[,1] <- trds
ytbl[names(frds),2] <- frds
omar <- par(mar=c(6,3.5,3,2))
barplot(t(ytbl),beside=TRUE,col=cls[1:2],ylim=c(0,max(ytbl)*1.2),las=2,
     names.arg=rownames(ytbl),cex.names=0.6,cex.axis=0.8)
abline(h=0)
legend("top",horiz=TRUE,fill=cls,legend=colnames(ytbl),cex=0.8)
par(mar=omar)
barplot(as.vector(ytbl[,2]/ytbl[,1]*100),col="lavender",ylim=c(0,100),
        names.arg=rownames(ytbl),las=2,cex.names=0.6,cex.axis=0.8)
abline(h=0)
title(main="Global yield",line=1)

dev.off()


# ###  Global yield by step

# load(file=file.path(repDir,"FLASH_table.RData"))
# rownames(flash.res) <- str_extract(rownames(flash.res),"^[A-Za-z0-9\\.-]+")

# p.nms <- rownames(PoolTbl)
# flash.res <- flash.res[p.nms,]
# rds.raw <- flash.res[,1]+flash.res[,2]
# rds.flash <- flash.res[,1]
# rds.MID <- PoolTbl[p.nms,"MIDReads"]
# rds.dmult <- PoolTbl[p.nms,"PrimerReads"]

# frdf.pool <- sapply(1:nrow(frdf), function(i)
               # samples$Pool.Nm[ which(samples$Patient.ID==frdf$Pat.ID[i] &
		                              # samples$Primer.ID==frdf$Ampl.Nm[i])[1] ])
								  
# rds.filt <- tapply(frdf$all,frdf.pool,sum)[p.nms]
# rds.ints <- tapply(frdf$Fn.rd,frdf.pool,sum)[p.nms]
# gbly <- cbind(raw=rds.raw,flash=rds.flash,MID=rds.MID,primer=rds.dmult,
              # filt=rds.filt,ints=rds.ints)
# Tgbl <- colSums(gbly)
# pct.gbl <- round(Tgbl[-1]/Tgbl[-6]*100,2)

# pdf.flnm3 <- "GlobalYieldBarplots.pdf"
# pdf(file.path(repDir,pdf.flnm3),paper="a4",width=6.5,height=10)
# par(mfrow=c(2,1))

# pal1 <- brewer.pal(8,"Dark2")
# pal2 <- brewer.pal(8,"Pastel2")
# m <- ncol(gbly)
# ymx <- max(gbly)*1.2
# bp <- barplot(t(gbly),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              # ylim=c(0,ymx),xaxt="n",cex.axis=0.8,ylab="# reads")
# axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
# abline(h=0)	
# legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
# title(main="Yield on pools by analysis step")

# bp <- barplot(t(gbly/gbly[,1]*100),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              # ylim=c(0,117),xaxt="n",cex.axis=0.8,ylab="Percentage")
# axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
# grid(nx=NA,ny=NULL)	
# abline(h=0)	
# legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
# title(main="Yield on pools by analysis step")

# par(mar=c(5,8,4,6))

# bp <- barplot(pct.gbl,col="lavender",border="navy",ylim=c(0,max(pct.gbl)),
        # ylab="yield (%)")
# text(bp,50,pct.gbl,col="navy",font=2,cex=0.8)
# title(main="Global yield by step")		

# boxplot(frdf$Fn.rd,border="gray",ylab="# of reads",outline=FALSE,
        # ylim=range(frdf$Fn.rd))
# points(jitter(rep(1,nrow(frdf)),a=0.10),frdf$Fn.rd,pch="+",cex=0.8)
# title(main="Final coverage")

# dev.off()
			  
# txt.flnm2 <- "GlobalYield-SumRprt.txt"
# sink(file=file.path(repDir,txt.flnm2))

# cat("\n   Global yield by analysis step")
# cat("\n===================================\n")
# cat("\nIn number of reads:\n\n")
# gbly <- rbind(gbly,TOTAL=colSums(gbly))
# print(gbly)	
		
# cat("\nIn percentage by step:\n\n")
# print(round(gbly[,-1]/gbly[,-6]*100,2))

# cat("\nIn percentage referred to raw reads:\n\n")
# print(round(gbly/gbly[,1]*100,2))			
# sink()
