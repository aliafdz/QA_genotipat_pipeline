
###################################################################
###  AQUESTA VERSIÓ RETALLA ELS PRIMERS PER AMBDOS EXTREMS 
###   ANANT-LOS A BUSCAR, SENSE PREFIXAR L'AMPLADA DE L'AMPLICÓ.
###################################################################

library(Biostrings)
library(ShortRead)

#-----------------------------------------------------------------------#

###  Formatar sencers emplenant amb 0 per l'esquerra
# Funció per ...
# 'nchar()' conta el nº de caràcters de l'argument
# 'substring()' extrau o reemplaça substrings d'un vector de caràcter
zeroFillInt2Char <- function(x,ln)
{ x <- paste("000000",x,sep="")
  substring(x,nchar(x)-ln+1)
}

# Funció per guardar la seqüència indicada en format fasta en la ruta de l'argument
write.fasta <- function(seqs,flnm)
{
  writeXStringSet(seqs,flnm) # Nota: es podria utilitzar aquesta funció directament
}

# Funció que dóna nom als haplotips en funció del nombre de mutacions que presenten respecte
# la seq de referència (que pot ser la màster) i de la seva freqüència poblacional, i guarda 
# les seqüències en un fitxer fasta.
SaveAllHaplotypes <- function(bseqs,nr,flnm)
{
  code <- "Hpl"
  ## Determinació de diferències respecte la seqüència màster (amb major freqüència)
  # 'pairwiseAlignment()' realitza un aliniament global de Needleman-Wunsch
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[1])
  # Recompte del nº de mismatch entre les seqüències
  nm <- nmismatch(psa)
  # Guarda les vegades que apareix cada nº de mismatches
  tnm <- table(nm)

  # Ordena els indexs de les seqüències per nombre de mutacions respecte la màster
  o <- order(nm)
  # Ordena les seqs segons el seu nº de mutacions
  bseqs <- bseqs[o]
  # Ordena també les freqüències segons el nº de mutacions de la seqüència
  nr <- nr[o]
  # Ordena la variable amb el nº de mismatches segons les mutacions (ordre ascendent)
  nm <- nm[o]

  ##  Numero d'ordre dins de cada nombre de mutacions
  # length(tnm) 
  isq <- unlist(sapply(1:length(tnm),function(i) 1:tnm[i]))
  
  ##  Ordenar per frequencia descendent dins de cada nombre de mutacions
  for(i in as.integer(names(tnm)))
  { idx <- which(nm==i)
    o <- order(nr[idx],decreasing=TRUE)
    bseqs[idx] <- bseqs[idx[o]]
    nr[idx] <- nr[idx[o]]
  }

  ##  Calcular frequencia relativa
  frq <- round(nr/sum(nr)*100,2)

  ## Nom complet per cada haplotipus
  nms <- paste(code,nm,zeroFillInt2Char(isq,4),sep=".")

  ##  Capçalera fasta amb nom, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ##  Salvar a fasta
  write.fasta(DNAStringSet(bseqs),file.path(trimDir,flnm))

  list(bseqs=bseqs,nr=nr,nm=nm)
}

#-----------------------------------------------------------------------#

###  Llegeix l'estructura de descripció de mostres
samples <- read.table(file.path(dataDir,"samples.csv"), sep="\t", header=T,
                      stringsAsFactors=F)
###  Llegeix el fitxer amb els primers
primers <- read.table(file.path(dataDir,"primers.csv"), sep="\t", header=T,
                    stringsAsFactors=F)

### Els primers específics emprats en l'amplificació (que afegeixen la cua M13 per seqüenciar) s'han d'eliminar
# dels amplicons a estudiar, ja que a l'emprar el mateix primer per totes les mostres es perden tots els possibles
# mismatches que determinen variabilitat genètica. 
# Per tant, es calculen les posicions de tall pels primers específics FW i RV. A la posició 5' del primer FW se li
# suma la seva longitud per determinar el punt de tall 5' de l'amplicó. 
# A la posició 5' del primer RV se li resta la longitud per determinar el punt de tall 3' de l'amplicó.
primers$FW.tpos <- primers$FW.pos+nchar(primers$Primer.FW)
primers$RV.tpos <- primers$RV.pos-nchar(primers$Primer.RV)

### Inicialitzacions
# Calcula el nº de mostres (files del fitxer samples)
Ns <- nrow(samples)
# Calcula el doble del total de mostres 
Ns <- Ns*2					 
# Genera una matriu amb el doble de files que el total de mostres i 4 columnes
pr.res <- matrix(0,nrow=Ns,ncol=4)
# Asigna el nom a les columnes de la matriu on s'aniran afegint els resultats
colnames(pr.res) <-  c("Tot.reads","matches","shorts","fn.reads")

# Genera un data frame amb el doble de files que el total de mostres i 9 columnes amb els noms indicats 
# 'character()' aplicat sobre un nombre enter retorna tants caràcters buits com el nº indiqui
# 'integer()' sobre un nombre enter retorna tants 0s com el nº indiqui
FlTbl <- data.frame(File.Name=character(Ns),Pat.ID=character(Ns),
                    Ampl.Nm=character(Ns),Pr.ID=integer(Ns),
					Str=character(Ns),Pos=integer(Ns),Len=integer(Ns),
					Reads=integer(Ns),Hpls=integer(Ns),
					stringsAsFactors=FALSE)
					

# Guarda el nom dels pools
pools <- unique(samples$Pool.Nm)
# Guarda un vector amb tants 0s com pools tenim
p.cv <- p.ok <- numeric(length(pools))
# Assigna al vector generat el nom dels pools
names(p.cv) <- names(p.ok) <- pools

# Genera els fitxers resultants en format .txt i .pdf que es guardaran a la carpeta reports
sink(file.path(repDir,"AmpliconLengthsRprt.txt"))
pdf(file.path(repDir,"AmpliconLengthsPlot.pdf"),paper="a4",
              width="6",height=10.5)
par(mfrow=c(2,1))

### Loop sobre el total de pools emprats (extrets del fitxer samples)
k <- 0
for(i in 1:length(pools))
{ # Identifica les mostres que corresponen al pool avaluat
  idx <- which(samples$Pool.Nm==pools[i])
  # Genera el nom dels fitxers .fna per cadascun dels MIDS que corresponen a les mostres indicades
  flnms <- paste(pools[i],".MID",samples$MID[idx],".fna",sep="")
  
  ### Loop sobre les mostres del pool avaluat -> cadascuna amb un MID concret
  for(j in 1:length(idx))
  { # Guarda l'index de la mostra avaluada
    jj <- idx[j]
  
    ### Carrega el fitxer fastq de la carpeta split corresponent al MID de la mostra avaluada
    # Si el fitxer no es troba al directori indicat, executa la següent iteració
    if(! file.exists(file.path(splitDir,flnms[j]))) next 
    # Llegeix el fitxer .fna del MID de la mostra avaluada
	  seqs <- readDNAStringSet(file.path(splitDir,flnms[j]))
	  # Suma el nº de reads del fitxer del MID al total del pool al qual correspon
	  p.cv[i] <- p.cv[i] + length(seqs)
    ### Retorna la posició en la columna Ampl.Nm del primer ID que s'està avaluant
	  # És a dir, de la mostra dins del pool que estem avaluant, indica quin és el seu primer ID i el busca
	  # al fitxer de primers
    ipr <- grep(samples$Primer.ID[jj],primers$Ampl.Nm)
  
    # Concatena els caràcters d'associació entre el pool que s'està avaluant, el MID i la regió de l'amplicó
	  tt <- paste("Pool",pools[i],"  MID",samples$MID[jj],
	            "  Ampl",primers$Ampl.Nm[ipr])
	  cat("\n\n",tt,"\n",sep="")
    		
    ### Elimina les seqüències del fitxer del MID avaluat que tinguin longitud menor a 180
	  seqs <- seqs[width(seqs)>180]

    ###  primer up matches 
    # Guarda la seq del primer FW específic de la regió avaluada
    pr.up <- primers$Primer.FW[ ipr ]
    # Busca la seq del primer FW en la regió 5' (posicions 1-100 definides al fitxer de paràmetres)
    # de les seqüències del fitxer .fna del MID avaluat. max.prdif defineix el màxim de mismatches permesos
    # Guarda totes les coincidències amb les posicions inicial, final, i longitud
    up.matches <- vmatchPattern(pattern=pr.up,
                            subject=subseq(seqs,start=target.io,end=target.in),
                            max.mismatch=max.prdif,fixed=FALSE)
    # Resta 1 a la posició inicial de cerca del primer (per defecte, 1)
    delta <- target.io-1
    # Aplica la funció que està definida al fitxer principal
    # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer FW 
    flags <- elementLengths(up.matches)>=1
    # Variable buida de nombre enter
    nr <- integer()
    # Variable buida de nombre enter
    shorts <- integer()
	  trim.len <- 0
	  # Genera un fitxer .fna el nom del qual està format pel ID del pacient, la regió (amplicó) avaluada,
	  # i la consecució PrFW.fna (primer forward)
    up.flnm <- paste(samples$Patient.ID[jj],primers$Ampl.Nm[ipr],"PrFW.fna",
	                 sep=".")
	  
    seqs.up <- ""
    if(sum(flags)) # Sumatori de seqüències amb 1 o més coincidències amb la seq del primer
    { # 'startIndex()' retorna una llista amb les posicions inicials de les coincidències del patró cercat
      # Cal tenir cura ja que si no hi ha coincidencia la llista guarda un valor NULL
      # Aplicar el condicional flags sobre la llista de posicions inicials permet eliminar els valors NULL
      # Guarda en una variable tots els valors de posició inicial del primer FW sobre les seqüències
      pos <- sapply(startIndex(up.matches)[flags],function(x) x[[1]])+delta
      
      ### Guarda les seqs amb coincidències del primer FW
      seqs.up <- seqs[flags]
	    ### Actualitza el nº de seqs amb les que no presentaven coincidències
	    seqs <- seqs[!flags]
	    
	    ### Retalla el primer 5' 
	    # Suma a la posició inicial del primer FW en la seq la seva longitud
      st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])
      # Retalla les seqüències des de la posició on acaba el primer FW fins al final 
      seqs.up <- subseq(seqs.up,start=st,end=width(seqs.up))   

      ### Localitza el primer per l'altre extrem
      # Guarda la reversa complementària del primer RV de la regió avaluada
      pr.p3 <- as.character(
	             reverseComplement(DNAString(primers$Primer.RV[ipr])))
      # Calcula la longitud de les seqüències que tenien coincidència amb el primer FW
      endp <- width(seqs.up)
      # Resta a la longitud de les seqs la posició final on hem de buscar el primer
      fstp <- endp-target.in
      # Calcula la diferència entre la nova posició final i la posició inicial
      delta <- fstp-1
      # Busca la seq del primer RV en la regió 3' (100 últimes posicions) de les seqüències on
      # ja s'ha retallat el primer FW. max.prdif defineix el màxim de mismatches permesos
      # Guarda totes les coincidències amb les posicions inicial, final, i longitud
      p3.matches <- vmatchPattern(pattern=pr.p3,
                            subject=subseq(seqs.up,start=fstp,end=endp),
                            max.mismatch=max.prdif,fixed=FALSE)
      # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer RV
      flags <- elementLengths(p3.matches)>=1
      
    ### Els reads que passen aquí tenen l'amplicó sencer, ja que han coincidit els dos primers
	  if(sum(flags))
	  { ### Retalla el primer RV per deixar l'amplicó net
	    seqs.up <- seqs.up[flags]
	    # Torna a fer el procés d'abans i guarda tots els valors de posició inicial del primer RV 
	    # sobre les seqüències
	    pos <- sapply(startIndex(p3.matches)[flags],function(x) x[[1]])+
	                  delta[flags] 
	    # Retalla les seqs des de la posició 1 fins on es detecta el primer RV
		  seqs.up <- subseq(seqs.up,start=1,end=pos-1)           
		
	    ### Colapsa les seqüències en haplotips + freqüències
		  # Ordena en ordre descendent les seqüències en funció de la seva freqüència (calculada amb 'table()')
      sqtbl <- sort(table(as.character(seqs.up)),decreasing=TRUE)
      # Guarda les seqüències en una variable
      bseqs <- names(sqtbl)
      # Els noms d'aquesta variable seran nombre de l'1 al total de seqs disponibles 
      names(bseqs) <- 1:length(bseqs)
      # Guarda les freqüències ordenades de les seqs
      nr <- as.integer(sqtbl)
      # Aplica la funció del principi sobre les seqüències, les seves freqüències i el fitxer .fna del 
      # pacient al que corresponen aquestes seqs (on es guardaran els resultats)
      lst <- SaveAllHaplotypes(bseqs,nr,up.flnm)
		  nr <- lst$nr
		
		cat("\nForward seqs, table of read lengths (over 10 rd)\n")
        tbl.len <- tapply(nr,nchar(bseqs),sum)
		print(tbl.len[tbl.len>=10])
		plot(as.integer(names(tbl.len)),tbl.len,type="h",
		     xlab="Read length",ylab="# reads")
        title(main=paste(tt," Str FW"))			 
      }
    }
	k <- k+1
	FlTbl$File.Name[k] <- up.flnm
	FlTbl$Pat.ID[k] <- samples$Patient.ID[jj]
	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
	FlTbl$Pr.ID[k] <- ipr
	FlTbl$Str[k] <- "fw"
	FlTbl$Pos[k] <- primers$FW.tpos[ipr]
	FlTbl$Len[k] <- mean(width(seqs.up))
    FlTbl$Reads[k] <- sum(nr)				   
    FlTbl$Hpls[k] <- length(nr)				   

	p.ok[i] <- p.ok[i] + sum(nr)

	pr.res[k,1] <- length(seqs)
    pr.res[k,2] <- sum(flags)
    pr.res[k,3] <- sum(shorts)
    pr.res[k,4] <- sum(nr)
  
    ###  primer dn matches 
    #################################
    pr.dn <- primers$Primer.RV[ ipr ]
    dn.matches <- vmatchPattern(pattern=pr.dn,
                            subject=subseq(seqs,start=target.io,end=target.in),
                            max.mismatch=max.prdif,fixed=FALSE)
    delta <- target.io-1
    flags <- elementLengths(dn.matches)>=1
    nr <- integer()
    shorts <- integer()
	
    dn.flnm <- paste(samples$Patient.ID[jj],primers$Ampl.Nm[ipr],"PrRV.fna",
	                 sep=".")
	seqs.dn <- ""
    if(sum(flags))
    { ###  Matches
      seqs.dn <- seqs[flags]
	  ###  Update remaining seqs
	  seqs <- seqs[!flags]
	  ###  Trim primer
	  pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta
      st <- pos + (primers$RV.pos[ipr]-primers$RV.tpos[ipr])   ###
      seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))   ###
	  seqs.dn <- seqs.dn[width(seqs.dn)>target.in+5]
	  
      ###  Reverse complement down matches
      seqs.dn <- reverseComplement(seqs.dn)
      ###  ... and match up primer
      io <- max(target.io-5,1)
      delta <- io-1
      dn.matches <- vmatchPattern(pattern=pr.up,
                              subject=subseq(seqs.dn,start=io,end=target.in+5),
                              max.mismatch=max.prdif,fixed=FALSE)
      flags <- elementLengths(dn.matches)>=1

      ###  Els que passen aquí tenen l'amplicó sencer
      if(sum(flags))
	  { ###  Trim reversed dn seqs
        seqs.dn <- seqs.dn[flags]
        pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta
        st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])   ###
        seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))   ###
		
		###  Collapse sequences to haplotypes+frequencies
        sqtbl <- sort(table(as.character(seqs.dn)),decreasing=TRUE)
        bseqs <- names(sqtbl)                
        names(bseqs) <- 1:length(bseqs)
        nr <- as.integer(sqtbl)
	    lst <- SaveAllHaplotypes(bseqs,nr,dn.flnm)
        nr <- lst$nr

        cat("\nReverse seqs, table of read lengths (over 10 rd)\n")
        tbl.len <- tapply(nr,nchar(bseqs),sum)
        print(tbl.len[tbl.len>=10])
		plot(as.integer(names(tbl.len)),tbl.len,type="h",
		     xlab="Read length",ylab="# reads")
        title(main=paste(tt," Str RV"))			 
	  }
	}  

    k <- k+1
	FlTbl$File.Name[k] <- dn.flnm
	FlTbl$Pat.ID[k] <- samples$Patient.ID[jj]
	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
	FlTbl$Pr.ID[k] <- ipr
	FlTbl$Str[k] <- "rv"
	FlTbl$Pos[k] <- primers$FW.tpos[ipr]
	FlTbl$Len[k] <- trim.len
    FlTbl$Reads[k] <- sum(nr)				   
    FlTbl$Hpls[k] <- length(nr)				   
	
	p.ok[i] <- p.ok[i] + sum(nr)

	pr.res[k,1] <- length(seqs)
    pr.res[k,2] <- sum(flags)
    pr.res[k,3] <- sum(shorts)
    pr.res[k,4] <- sum(nr)
  } # Fi del loop sobre totes les mostres corresponents al pool avaluat
} # Fi del loop sobre els dos pools avaluats
sink()
dev.off()

###  Sincronitzar structures
fl <- FlTbl$Pr.ID>0
FlTbl <- FlTbl[fl,]
pr.res <- pr.res[fl,]
anms <- paste(FlTbl$Pat.ID,FlTbl$Ampl.Nm,sep=".")

rownames(samples) <- paste(samples$Patient.ID,samples$Primer.ID,sep=".")

###  Plot results
library(RColorBrewer)
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")

fw.idx <- which(FlTbl$Str=="fw")
rv.idx <- which(FlTbl$Str=="rv")
mprres <- data.frame(PatID=FlTbl$Pat.ID[fw.idx],
                     PrimerID=FlTbl$Ampl.Nm[fw.idx],
                     Treads=pr.res[fw.idx,1],
					 Shorts=pr.res[fw.idx,3]+pr.res[rv.idx,3],
                     FW.match=pr.res[fw.idx,4],
                     RV.match=pr.res[rv.idx,4],
					 Fn.reads=pr.res[fw.idx,4]+pr.res[rv.idx,4],
                     stringsAsFactors=FALSE)
mres <- mprres
T.reads <- apply(mres[,5:6],2, function(x)
              tapply(x,mres$PatID,sum))
if(length(unique(mres$PatID))==1)
{ x <- matrix(T.reads,nrow=1)
  rownames(x) <- mres$PatID[1]
  colnames(x) <- names(T.reads)
  T.reads <- x
}  
  
pdf.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=11)
par(mfrow=c(2,1),mar=c(7,4,4,2)+0.1)

ymx <- max(rowSums(T.reads))*1.2
barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

res.mat <- mres[,5:6]
ymx <- max(res.mat)*1.2
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
dev.off()

yield <- p.ok/p.cv*100
PoolTbl <- data.frame(MIDReads=p.cv,PrimerReads=p.ok,Pct=yield)

txt.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.txt",sep="_")
sink(file.path(repDir,txt.flnm))
cat("\nTable of reads identified by primer\n\n")
print(mprres)
cat("\n")
print(FlTbl[,-1])
cat("\nTotal reads identified by patient\n\n")
print(T.reads)
cat("\nYield by pool\n\n")
print(PoolTbl)
sink()

###  Save tables

# FlTbl <- FlTbl[FlTbl$Reads>0, ]  # Erase rows of null files
save(FlTbl,PoolTbl,file=file.path(repDir,"SplittedReadsFileTable.RData"))
sink(file.path(repDir,"SplittedReadsFileTable.txt"))
print(FlTbl)
sink()


###  Plot on A4 horizontally
pdf.flnm <- paste(proj.nm,"SplitByPrimersOnFlash-hz.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4r",width=10,height=6)
par(mar=c(7,4,4,2)+0.1)

ymx <- max(res.mat)*1.2
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))

ymx <- max(rowSums(T.reads))*1.2
bp <- barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx),xaxt="n")
axis(side=1,at=bp,rownames(T.reads),cex.axis=0.6,las=2)
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

par(mfrow=c(1,2))

ymx <- max(c(res.mat[,1],res.mat[,2]))
rownames(primers) <- primers$Ampl.Nm
reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="fw"],"Region"]
boxplot(res.mat[,1]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,1],
       pch="+",cex=0.8)
title(main="Primers identified on forward reads")	   

reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="rv"],"Region"]
boxplot(res.mat[,2]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,2],
       pch="+",cex=0.8)
title(main="Primers identified on reverse reads")	   

dev.off()
