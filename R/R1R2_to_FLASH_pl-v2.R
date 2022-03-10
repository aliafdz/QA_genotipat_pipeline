# Defineix una funció que s'aplicarà més endavant
fastq.nreads <- function(flnm)
{
  n <- 0
  ###  Aplica streamer (iteració) en el fitxer fastq 
  strm <- FastqStreamer(flnm,n=chunck.sz) # chunck.sz definit en el fitxer principal
  ###  Carrega el fitxer fastq per chuncks i actualitza en nº total de reads
  while(length(sqq <- yield(strm)))
    n <- n + length(sqq)	
  close(strm)
  return(n)
}

#--------------------------------------------------------------------#

# Llista de fitxers disponibles en la carpeta run
flnms <- list.files(runDir)
# La funció sub() permet substituir un patró pel que indiquem com 2n argument
# En aquest cas, les variables flnms i snms son idèntiques (de moment)
snms <- sub("\\.fastq$","",flnms)

# Taula on es disposa els nom de cada fitxer fastq en la primera columna i després es 
# separa cadascun del termes del nom del fitxer en successives columnes
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
# Elimina les columnes 3 i 5 de la taula (que indiquen L001 i l'extensió)
parts <- parts[,-c(3,5)]
# Assigna els noms a les columnes: ID pacient, ID mostra/pool, read
colnames(parts) <- c("PatID","SmplID","Read")

# Guarda a la variable la ruta dels fitxers que es troben a la carpeta run
flnms <- file.path(runDir,flnms)
# Dels 4 fitxers que hi havia, guarda en 2 variables els que corresponen a R1 i R2
R1.flnms <- flnms[parts[,3]=="R1"] # parts[,3] = columna 3 (read) de la taula parts
R2.flnms <- flnms[parts[,3]=="R2"]

# Guarda dos noms de fitxer corresponents a R1: pool, ID de mostra i extensió flash.fastq
# Com hi ha 2 fitxers amb R1, vol dir que avaluem 2 regions o pools
out.flnms <- paste(parts[parts[,3]=="R1",1],parts[parts[,3]=="R1",2],
                   "flash.fastq",sep="_")
# Genera la ruta dels fitxers definits abans a la carpeta flash
out.flnms <- file.path(flashDir,out.flnms)

# Guarda de la taula només els R1
parts <- parts[parts[,3]=="R1",,drop=FALSE]
# Construeix una matriu 2x2 amb tot 0
res <- matrix(0,nrow=length(out.flnms),ncol=2) # length(out.flnms)= 2 en aquest cas, que son els pools
# Assigna com a nom de fila el pool (regió de VHB) i com a columnes els reads segons si s'ha donat o no extensió
rownames(res) <- paste(parts[,1],parts[,2],sep="_")
colnames(res) <- c("Extended","NoExtd")

# Itera sobre el nº de fitxers R1 guardats abans (nº de pools)
for(i in 1:length(out.flnms))
{ 
  # Comprova que estiguin els dos fitxers R1 i R2 de la iteració (pool) avaluada
  ok <- file.exists(R1.flnms[i]) & file.exists(R2.flnms[i])
  if(!ok) next # Si no es compleix la condició continua la següent iteració
  # Concatena la ruta de l'executable flash, els paràmetres definits en el fitxer QA i la ruta dels fitxers
  # R1 i R2 del pool avaluat
  command <- paste(flash,flash.opts,R1.flnms[i],R2.flnms[i],collapse=" ")
  # La funció system() invoca una comanda especificada en l'argument
  # Executa el programa FLASH -> es guarden els fitxers següents a la carpeta global: out.extendedFrags.fastq, 
  # out.hist, out.histogram, out.notCombined_1.fastq i out.notCombined_2.fastq
  es <- system(command,intern=FALSE,wait=TRUE,show.output.on.console=FALSE,
              ignore.stdout=FALSE,invisible=TRUE)
  if(es!=0) next
  # Copia el fitxer de 'from' en 'to'. Es guarda el fitxer resultant (to) en la carpeta flash
  file.copy(from="out.extendedFrags.fastq",to=out.flnms[i],overwrite=TRUE)
  
  # sqq <- readFastq(dirPath=".",patt="out.extendedFrags.fastq")
  # res[i,1] <- length(sqq)
  ## Llegeix el fitxer out.hist, una taula que recull els valors de 2 variables V1 i V2
  hstln <- read.table("out.hist",header=FALSE)
  # De la matriu res d'abans, que estava buida, recull en la columna 1 (extended) + fila del pool avaluat
  # el sumatori de la V2 de la variable anterior
  res[i,1] <- sum(hstln[,2])

  #sqq <- readFastq(dirPath=".",patt="out.notCombined_1.fastq")
  #res[i,2] <- length(sqq)
  ## Per la 2a columna de la matriu (no extended), aplica la funció definida al principi, que retorna el 
  # total de reads que hi ha en el fitxer de l'argument
  res[i,2] <- fastq.nreads("out.notCombined_1.fastq")
}

# Guarda un dataframe amb les dades resultants del FLASH (taula res), afegint la 
# columna Yield calculada dividint els reads extended entre el total *100
df.res <- data.frame(res,Yield=round(res[,1]/(res[,1]+res[,2])*100,1))
# Guarda el fitxer de report de FLASH en format .txt
txt.flnm <- file.path(repDir,"FLASH_report.txt")
# Comandes per omplir el fitxer txt generat. Recull els paràmetres indicats al FLASH
# i la taula df.res que conté els resultats obtinguts
sink(txt.flnm)
cat("\nExtending Illumina reads by FLASH\n")
cat("\nFLASH parameters:")
cat("\n    Minimum overlap:",min.ov)
cat("\n    Maximum overlap:",max.ov)
cat("\n        Error level:",err.lv,"\n\n")
print(df.res)
sink() # Tanca el fitxer

flash.res <- df.res # Necessari??
# Guarda a la carpeta reports la taula en format RData
save(flash.res,file=file.path(repDir,"FLASH_table.RData"))

# Genera el pdf que contindrà el gràfic barplot dels resultats FLASH
pdf.flnm <- file.path(repDir,"FLASH_barplot.pdf")
pdf(pdf.flnm,paper="a4",width=6,height=10)
par(mfrow=c(2,1))

library(RColorBrewer)

pal=brewer.pal(8,"Dark2") # Crea la paleta de colors
M <- data.matrix(df.res[,1:2]) # Guarda les columnes 1 i 2 de la taula de resultats
ymx <- max(M)*1.15 # Defineix límit superior del gràfic
barplot(t(M),beside=TRUE,las=2,col=pal[1:2],ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal[1:2],cex=0.8,
       legend=c("Extended","Not extended"))
title(main="FLASH results on paired-ends",line=2.5)
title(main="Yield in number of reads",cex.main=1,line=1)

# Genera un altre barplot on s'inclouen els resultats del Yield de FLASH
bp <- barplot(df.res$Yield,col="Lavender",ylim=c(0,100),ylab="Yield",
              names.arg=rownames(df.res),las=2)	   
text(bp,10,paste(df.res$Yield,"%",sep=""),srt=90,col="navy")
title(main="Yield in percentage",cex.main=1,line=1)

dev.off()
