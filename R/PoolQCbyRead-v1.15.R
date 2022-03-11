
library(ShortRead)
library(Biostrings)
library(stringr)

# Funció amb el mateix nom que en el fitxer PoolQCbyPos però acaba diferent
fn.fastq <- function(flnm,ln=301)
{ # Variable numèrica
  nrds <- numeric()
  # Variable de nombre enter
  all.ln <- integer()
  # Variable de nombre enter
  all.nl30 <- integer()
  
  ### Aplica streamer (iteració) en el fitxer fastq
  strm <- FastqStreamer(flnm,n=chunck.sz) # chunck.sz definit en el fitxer principal
  
  ### Carrega el fitxer fastq per chuncks
  nchk <- 0
  while(length(sqq <- yield(strm)))
  { nchk <- nchk+1 # Nº de cada iteració
    nrds[nchk] <- length(sqq) # L'item de la iteració avaluada correspondrà al nº de reads
    
    ### Phred scores. Codi ASCII-33
    # Funció 'quality()' retorna el valor de qualitat dels strings
    # Funció 'as()' permet fer coerció del resultat a matriu
    phrsc <- as(quality(sqq),"matrix")
    
    ### A la matriu amb els Phred scores aplica el sumatori de les bases amb puntuació menor a 30 (per sota de Q30)
    nl30 <- apply(phrsc,1,function(x) sum(x<30,na.rm=TRUE)) # na.rm per no tenir en compte missing values
    # La variable de nombre enter ara inclourà les bases per sota de Q30
    all.nl30 <- c(all.nl30,nl30)
    ###  Longituds de seqüència
    sqln <- width(sqq) # Cicles de seqüenciació 	
    all.ln <- c(all.ln,sqln) # Vector amb el nº de cicles com a nombres enters
  }
  close(strm)
  # Retorna una llista formada per les 2 matrius: amb el nº de cicles i les bases per sota de Q30
  return(list(all.ln=all.ln,all.nl30=all.nl30))
}  

#-------------------------------------------------------------------------#

flashDir <- "./flash"
repDir <- "./reports"
dataDir <- "./data"

### Fitxers que resulten de Flash				
flnms <- list.files(flashDir)
# Guarda del fitxer només el nom del pool seguit de S1 o S2
snms <- sub("_flash\\.fastq$","",flnms)

# Genera una taula amb el nom del pool en una columna i S1 o S2 en una altra
parts <- t(sapply(snms,function(str) strsplit(str,split="_")[[1]]))
if(is.vector(parts))
  parts <- matrix(parts,nrow=1)
colnames(parts) <- c("PatID","SmplID")

### Loop sobre pools
for(i in 1:length(snms))
{ # Aplica la funció del principi sobre el fitxer dins de flash que correspongui
  lst1 <- fn.fastq(file.path(flashDir,flnms[i]))
  # Genera el pdf on aniran els gràfics
  pdf.flnm <- paste("PoolQCbyRead_",parts[i,"PatID"],".pdf",sep="")
  pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=10.5)
  par(mfrow=c(2,1))
  # Després d'aplicar la funció, de la llista resultant agafa la matriu de les bases<Q30, i agafa només
  # aquelles per sota del quantil a 0.99
  nl30 <- lst1$all.nl30[lst1$all.nl30<quantile(lst1$all.nl30,p=0.99)]
 
  # Histograma amb el subset anterior
  hdt <- hist(nl30,breaks=30,main="",xlab="# bases below Q30")
  # Aplica els 3 quantils indicats a les bases<Q30
  qs <- quantile(lst1$all.nl30,p=c(0.5,0.8,0.95))
  # Dibuixa una fletxa blava en l'eix X per indicar el rang del quantil 0.95
  arrows(x0=0,x1=qs[3],y0=0,code=3,angle=90,len=0.1,col="blue",lwd=2)
  # Dibuixa una fletxa rosa en l'eix X per indicar el rang del quantil 0.80
  arrows(x0=0,x1=qs[2],y0=0,code=3,angle=90,len=0.1,col="maroon",lwd=2)
  # Dibuixa una línia vertical blava indicant el quantil 0.5
  abline(v=qs[1],lty=4,col="blue",lwd=2)
  q95 <- paste("q95% ",qs[3])
  q80 <- paste("q80% ",qs[2])
  med <- paste("q50% ",qs[1])
  # hist$mids= the n cell midpoints; hist$counts= n integers; for each cell, the number of x[] inside.
  text(x=max(hdt$mids)*0.95,y=max(hdt$counts)*0.9,adj=1,cex=0.8,
       paste(q95,q80,med,sep="\n"))
  title(main=parts[i,"PatID"],line=1)
  
  # Dividez les bases<30 entre el nº de cicles -> fracció de les bases per read
  all.fnl30 <- lst1$all.nl30/lst1$all.ln
  # D'aquesta matriu de divisió agafa els valors per sota del quantil 0.99
  fnl30 <- all.fnl30[all.fnl30<quantile(all.fnl30,p=0.99)]
  # Nou histograma però amb les fraccions de bases per read
  hdt <- hist(fnl30,breaks=30,main="",
              xlab="fraction of bases below Q30 by read")
  qs <- round(quantile(all.fnl30,p=c(0.5,0.8,0.95)),3)
  arrows(x0=0,x1=qs[3],y0=0,code=3,angle=90,len=0.1,col="blue",lwd=2)
  arrows(x0=0,x1=qs[2],y0=0,code=3,angle=90,len=0.1,col="maroon",lwd=2)
  abline(v=qs[1],lty=4,col="blue",lwd=2)
  q95 <- paste("q95%",qs[3])
  q80 <- paste("q80%",qs[2])
  med <- paste("q50%",qs[1])
  # Valors amb 3 decimals
  vals <- sprintf("%5.3f",c(qs[3],qs[2],qs[1]))
  vals <- paste(c("q95%","q80%","q50%"),vals,collapse="\n")
  text(x=max(hdt$mids)*0.98,y=max(hdt$counts)*0.95,adj=1,cex=0.8,vals)
  # Línies verticals grises en les posicions x indicades
  abline(v=c(0.01,0.02,0.05),lty=4,col="gray")
  
  # Calcula el sumatori (en %) de reads que disposen de menys de 1, 2 o 5% de les seves bases per sota de Q30
  # Si Pr(f<=0.05)=82%, per exemple, vol dir que si eliminem els reads amb més del 5% de les bases<Q30 encara
  # ens quedem amb un 82% dels reads totals. 
  p1pct <- sum(all.fnl30<=0.01)/length(all.fnl30)*100
  p2pct <- sum(all.fnl30<=0.02)/length(all.fnl30)*100
  p5pct <- sum(all.fnl30<=0.05)/length(all.fnl30)*100
  vals <- sprintf("%5.1f",c(p1pct,p2pct,p5pct))
  txt <- paste("Pr(f<=",c(1,2,5),"%)  ",vals,"%",sep="")
  txt <- paste(txt,collapse="\n")
  text(x=max(hdt$mids)*0.98,y=max(hdt$counts)*0.75,adj=1,cex=0.8,txt)
  title(main=parts[i,"PatID"],line=1)
  dev.off()
}
