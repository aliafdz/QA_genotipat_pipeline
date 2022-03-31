
#######################################################################
###
###   Demultiplexat de fitxers fastq molt grans per streamer
###
###   Simplement demultiplexa, no retalla res (els MIDs s'hi queden)
###   També descarta les seqs més curtes de 150
###
#######################################################################

library(Biostrings)
library(ShortRead)

max.middif <- 1

###  Llegeix l'estructura de descripció de mostres
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      stringsAsFactors=F)
# Guarda la columna Pools (regions 5' o preS1). 'unique' per eliminar elements duplicats (només hi haurà 2)
pools <- unique(samples$Pool.Nm)

###  Llegeix el fitxer d'identificació dels mids
mids <- read.table(file.path(dataDir,"mids.csv"), sep="\t", header=T,
                   stringsAsFactors=F)

###  Llista de fitxers fastq disponibles en la carpeta flashFilt (directori redefinit) 
flnms <- list.files(flashDir)

###  Inicialitzacions
# 'paste' per concatenar strings separats per punt.
# Així s'identifica cada MID (mostra) amb la regió avaluada (pool).
nms <- unique(paste(samples$Pool.Nm,samples$MID,sep="."))

# 'strsplit' separa la concatenació del pas anterior; 'as.data.frame' ho transforma en taula, però està en horitzontal.
# 't' fa la transposició perquè estigui en vertical (en files)
items <- t(as.data.frame(strsplit(nms,split="\\.")))

# Guarda les entrades de nms, serien el total de mostres (2 per pacient)
Ns <- length(nms)

# Taula d'abans (convertint els MIDS com a nombres enters) afegint una 3a columna anomenada Reads, que inicialment 
# és tot 0.
pr.res <- data.frame(Pool.Nm=items[,1],MID=as.integer(items[,2]),
                     Reads=integer(Ns),stringsAsFactors=FALSE)

# Indica com a nom de les files de la taula l'indicador nms (del tipus Pool.MID) 
rownames(pr.res) <- nms
# Guarda una llista de tants 0 com pools tenim, en aquest cas 2. Aquí es guardarà el nº de reads no assignats
MID0.reads <- integer(length(pools))
# Associa aquests 0 als 2 pools (li assigna el nom)
names(MID0.reads) <- pools
# Guarda la llista on es guardarà el nº de reads totals per pool
pool.reads <- integer(length(pools))
names(pool.reads) <- pools


###  Loop sobre pools en samples
for(i in 1:length(pools))
{
  ###  Identificar fastq del pool
  # Busca dins dels arxius de la carpeta flashDir el fitxer per cada pool (regió VHB).
  # 'grep' serveix per buscar un patró dins d'un vector de caracters, i retorna l'index del vector que encaixa amb la cerca
  ip <- grep(paste(pools[i],"_",sep=""),flnms) 
  # Si no troba l'arxiu de la iteració que està fent se'l salta i passa a la següent iteració.
  if(length(ip)==0) next
  
  ###  Identificar MIDs en pool
  # 'which' retorna els elements que compleixen el que s'indica com argument. En aquest cas, 
  # quins elements de la columna Pool coincideixen amb el que s'està iterant.
  idx <- which(samples$Pool.Nm==pools[i])
  # Igual que abans, si no coincideix cap torna a iterar. 
  if(length(idx)==0) next
  
  ###  Identificar MIDs en pool
  # Guarda els MIDS que corresponen al pool que s'està avaluant. 
  mid.set <- unique(samples$MID[idx])
  
  ###  Aplicar streamer en el fitxer fastq 
  # La funció 'FastqStreamer' llegeix l'arxiu fastq que estem iterant del directori flash i retorna
  # parts successives d'aquest.
  strm <- FastqStreamer(file.path(flashDir,flnms[ip]))
  p.cv <- 0 # nº inicial de reads a 0
  app.flag <- m0.app.flag <- FALSE
  
  ###  Carrega el fitxer fastq per chuncks
  # Itera sobre tots els blocs del fastq que inclouen: id, seq, +, qualitats
  while(length(sqq <- yield(strm)))
  { seqs <- sread(sqq) # Guarda la seq de cada conjunt o chunck
  
  ###  Actualitza el nº total de reads
  p.cv <- p.cv + length(seqs)
  ###  Es descarten aquells reads de longitud menor a 150
  seqs <- seqs[ width(seqs) > 150 ]
  
  ###  Loop sobre els MIDs del pool avaluat
  for(j in mid.set)
    # Guarda en una variable la seq inclosa en les posicions indicades.
    # En aquest cas, es busca el MID entre les posicions 1-40 (start-end), definides 
    # a l'arxiu de paràmetres. 
  { sbsq <- subseq(seqs,start=mid.start,end=mid.end)
  k <- which(mids$MID.ID==j) # Guarda els MIDs que coincideixen amb l'avaluat.
  pr.up <- mids$MID.Seq[k] # Guarda la seq del MID que ha coincidit. 
  # Recompte de coincidències entre la seq 1-40 extreta i la del MID que ha 
  # coincidit abans, amb màxim 1 diferència de mismacth.
  up.matches <- vcountPattern(pattern=pr.up,subject=sbsq,
                              max.mismatch=max.middif,fixed=TRUE) 
  # Només guarda els que tinguin més d'una coincidència. 
  flags <- up.matches>=1
  if(sum(flags)) # La funció 'sum' assegura que tinguem un valor major a 1.
    
    # Guarda en un arxiu .fna la seq dels reads que s'han associat al MID
    # que s'està avaluant, amb nom `pool.nºMID.fna`, en el directori splits.
  { KK <- which(nms==paste(pools[i],j,sep="."))[1]
  pr.res$Reads[KK] <-  pr.res$Reads[KK]+sum(flags)
  up.flnm <- paste(pools[i],".MID",j,".fna",
                   sep="")
  writeXStringSet(seqs[flags],file.path(splitDir,up.flnm),
                  append=app.flag)
  seqs <- seqs[!flags]  # Actualitza els reads que no han fet match en aquesta iteració.					  
  }
  }
  pool.reads[i] <- p.cv # Sumatori dels reads obtinguts en el pool avaluat.
  
  if(length(seqs)) # En cas que s'hagin iterat tots els MIDs i quedin reads sense assignar:
    # Guarda un altre arxiu .fna amb els reads que no s'han assignat a MID i ho anomena MID0.
  { mid0.flnm <- paste(pools[i],"MID0.fna",sep=".")
  writeXStringSet(seqs,file.path(splitDir,mid0.flnm),
                  append=m0.app.flag)
  m0.app.flag <- TRUE
  MID0.reads[i] <- MID0.reads[i]+length(seqs)
  }
  app.flag <- TRUE
  }
  close(strm)
}

# Guarda un data frame que inclogui la quantificació de reads totals, els que no s'han assignat (MID0)
# i els que sí s'han assignat
by.pools <- data.frame(TotReads=pool.reads,NoMID=MID0.reads,
                       MIDReads=tapply(pr.res$Reads,
                                       factor(pr.res$Pool.Nm,levels=pools),sum))
# Genera el pdf indicar que es guardarà a la carpeta reports
pdf.flnm <- file.path(repDir,"SplidByMIDs.barplots.pdf")
pdf(pdf.flnm,paper="a4",width=5.5,height=10)
par(mfrow=c(2,1),mar=c(7.5,4,4,2)+0.1)

library(RColorBrewer)
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
# Genera un gràfic de barres amb les dades del nº de reads per MID
bp <- barplot(pr.res$Reads,col="lavender",border="navy")
axis(1,at=bp,las=2,rownames(pr.res),cex.axis=0.8)
title(main="Coverage by Pool and MID",line=1.5)

# Defineix el límit de l'eix y a partir del màxim de reads totals obtinguts
ymx <- max(data.matrix(by.pools))*1.2
# Genera un altre gràfic de barres amb les dades del nº de reads per pool (assignats a MID i no assignats)
bp <- barplot(t(data.matrix(by.pools)),beside=TRUE,ylim=c(0,ymx),
              col=pal2[1:3],border=pal1[1:3])
legend("top",horiz=TRUE,fill=pal2[1:3],legend=colnames(by.pools),cex=0.8)			  
title(main="Coverage by Pool",line=1.5)

dev.off()
# Genera un fitxer .txt que es guardarà a la carpeta reports, on es guadarà en format taula els resultats
# representats als gràfics
txt.flnm <- file.path(repDir,"SplidByMIDs.Rprt.txt")
sink(txt.flnm)
cat("\nCoverage by Pool and MIDs\n\n")
rownames(pr.res) <- NULL
print(pr.res)
cat("\n")
print(by.pools)
sink()
