
###################################################################
###  AQUESTA VERSIÓ RETALLA ELS PRIMERS PER AMBDOS EXTREMS 
###   ANANT-LOS A BUSCAR, SENSE PREFIXAR L'AMPLADA DE L'AMPLICÓ.
###################################################################

library(Biostrings)
library(ShortRead)

#-----------------------------------------------------------------------#

# Funció per assignar a cada haplotip el seu nom en funció del nº de mutacions que presenta
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

# Funció per guardar la seqüència indicada en format fasta en la ruta de l'argument
write.fasta <- function(seqs,flnm)
{
  writeXStringSet(seqs,flnm) # Nota: es podria utilitzar aquesta funció directament
}

# Funció que dóna nom als haplotips en funció del nombre de mutacions que presenten respecte
# la seq de referència (que pot ser la màster) i de la seva freqüència poblacional, i guarda 
# les seqüències en un fitxer fasta.
SaveAllHaplotypes <- function(bseqs,nr,flnm)
{ # Aquesta variable determina el nom dels fitxers resultants per cada haplotip
  code <- "Hpl"
  ## Determinació de diferències respecte la seqüència màster (amb major freqüència)
  # 'pairwiseAlignment()' realitza un aliniament global de Needleman-Wunsch
  psa <- pairwiseAlignment(pattern=bseqs,subject=bseqs[1])
  # Recompte del nº de mismatch entre les seqüències
  nm <- nmismatch(psa)
  # Guarda les vegades que apareix cada nº de mismatches
  tnm <- table(nm)

  # Ordena el nombre de mutacions respecte la màster en ordre ascendent
  o <- order(nm)
  # Ordena les seqs segons el seu nº de mutacions respecte la màster
  bseqs <- bseqs[o]
  # Ordena també les freqüències segons el nº de mutacions de la seqüència
  nr <- nr[o]
  # Ordena la variable amb el nº de mismatches segons les mutacions (ordre ascendent)
  nm <- nm[o]

  ## Assigna el número d'ordre dins de cada nombre de mutacions:
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

  ## Calcula la freqüència relativa de cada amplicó pertanyent al MID avaluat
  frq <- round(nr/sum(nr)*100,2)

  ## Nom complet per cada haplotip
  # Defineix el nom que rebrà cadascun dels haplotips associats al MID avaluat
  # code està definit al principi (Hpl), nm= mutacions respecte seq màster
  # 'zeroFillInt2Char' definida al principi, retorna el nº d'ordenació de l'haplotip dins del conjunt
  # amb el mateix nº de mutacions
  # Exemple: "Hpl.2.0387" indica que es tracta de l'haplotip nº387 que presenta 2 mutacions
  nms <- paste(code,nm,zeroFillInt2Char(isq,4),sep=".")

  ## Capçalera del fitxer fasta amb nom de l'haplotip, nombre de reads i freq relativa
  names(bseqs) <- paste(nms,nr,frq,sep="|")

  ## Generació del fitxer fasta amb tots els haplotips assignats a la carpeta trim
  write.fasta(DNAStringSet(bseqs),file.path(trimDir,flnm))
  # Guarda una llista amb les seqüències, les seves freqüències i mutacions respecta la màster (en aquest ordre)
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

# Genera un data frame amb el doble de files que el total de mostres i 9 columnes amb els noms dels arguments 
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
  # Guarda el nom dels fitxers .fna dels MIDS corresponents a les mostres indicades, dels quals parteix
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
	  
	  ### Per guardar el nº total de reads abans de demultiplexar per primers:
	  # totalseqs <- length(seqs)
	  
	  # Afegeix el total de reads del fitxer del MID concret al total del pool al qual correspon
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

 ### Coincidències cadena forward
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
	  # Genera un fitxer .fna, el nom del qual està format per l'ID del pacient, la regió (amplicó) avaluada,
	  # i la consecució PrFW.fna (primer forward), que es guarda a la carpeta trim
    up.flnm <- paste(samples$Patient.ID[jj],primers$Ampl.Nm[ipr],"PrFW.fna",
	                 sep=".")
	  
    seqs.up <- "" # Coincidències cadena FW
    if(sum(flags)) # Sumatori de seqüències amb 1 o més coincidències amb la seq del primer
    { # 'startIndex()' retorna una llista amb les posicions inicials de les coincidències del patró cercat
      # Cal tenir cura, ja que si no hi ha coincidència la llista guarda un valor NULL
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

      ### Localitza el primer per l'altre extrem --> a la mateixa cadena! 
      # Guarda la reversa complementària del primer RV de la regió avaluada
      pr.p3 <- as.character(
	             reverseComplement(DNAString(primers$Primer.RV[ipr])))
      # Calcula la longitud de les seqüències que tenien coincidència amb el primer FW després de retallar-lo
        # Seria equivalent a dir que es calcula la posició final del primer RV (parlant en sentit 5'-3')
      endp <- width(seqs.up)
      # Resta a la longitud de les seqs la posició final on hem de buscar el primer
        # Com estem parlant del primer RV, que es troba a 3', l'haurem de buscar en les 100 últimes posicions
      fstp <- endp-target.in
      # Calcula la diferència entre la posició on comença el primer RV i la posició inicial de la seqüència
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
		
	    ### Colapsa les seqüències 'up' en haplotips + freqüències
		  # Ordena en ordre descendent les freqüències (calculades amb 'table()') de cada haplotip
		    # Les seqüències retallades es transformen en caràcters, i a l'aplicar la funció 'table()', es comparen
		    # entre elles, per tant totes les seqs idèntiques s'agrupen juntes. El nº de seqs idèntiques defineix la freq
      sqtbl <- sort(table(as.character(seqs.up)),decreasing=TRUE)
      # Guarda les seqüències en una variable
      bseqs <- names(sqtbl)
      # Els noms d'aquesta variable seran nombres de l'1 al total de seqs disponibles 
      names(bseqs) <- 1:length(bseqs)
      # Guarda les freqüències dels haplotips (com a nombres enters) ordenades de les seqs
      nr <- as.integer(sqtbl)
      # Aplica la funció del principi sobre les seqüències, les seves freqüències i el fitxer .fna del 
      # pacient al que corresponen aquestes seqs (on es guardaran els resultats)
      lst <- SaveAllHaplotypes(bseqs,nr,up.flnm)
      # Guarda les freqüències dels haplotips obtinguts
		  nr <- lst$nr
		
		
		cat("\nForward seqs, table of read lengths (over 10 rd)\n")
    # 'tapply()' en aquest cas calcula el sumatori (sum) de les freqüències (nº reads) segons les diferents
		# longituds de seqüència
		tbl.len <- tapply(nr,nchar(bseqs),sum)
		# Filtra les longituds amb més de 10 reads (per això el títol amb la funció 'cat()') 
		print(tbl.len[tbl.len>=10])
		# Gràfic representant les diferents longituds en X i les freqüències (nº reads) en Y
		plot(as.integer(names(tbl.len)),tbl.len,type="h", # "h" per representar histogrames en forma de línies verticals
		     xlab="Read length",ylab="# reads")
        title(main=paste(tt," Str FW"))			 
      }
    }
    # Actualitza k per indicar l'index a la taula de reports: per cada mostra (2 per pacient) hi haurà 2 resultats
    # corresponents a les cadenes forward i reverse
  	k <- k+1
  	# Nom del fitxer (pacient, regió, PrFw)
  	FlTbl$File.Name[k] <- up.flnm
  	# Identificació del pacient
  	FlTbl$Pat.ID[k] <- samples$Patient.ID[jj]
  	# Coordenades de la egió amplificada (5'X o preS1)
  	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
  	# Identificador del primer (1 o 2) segons la regió que amplifica
  	FlTbl$Pr.ID[k] <- ipr
  	# Cadena (en aquest cas forward)
  	FlTbl$Str[k] <- "fw"
  	# Posició del genoma de HBV en la que comença l'amplicó (després d'eliminar els primers)
  	FlTbl$Pos[k] <- primers$FW.tpos[ipr]
  	# Mitjana de longitud de les seqüències després de retallar-les
  	  # Comprovant els valors de width(seqs.up) es pot veure que no totes les seqüències tenen
  	  # la mateixa longitud, per tant es fa la mitjana de totes
  	FlTbl$Len[k] <- mean(width(seqs.up)) #### 
  	# Sumatori del total de reads en el MID avaluat (mostra) que s'han assignat
    FlTbl$Reads[k] <- sum(nr)				   
    # Total d'haplotips detectats
    FlTbl$Hpls[k] <- length(nr)				   
    
    # Afegeix el nº de reads (del MID concret) idenficats amb la cadena FW al total del pool al qual correspon
  	p.ok[i] <- p.ok[i] + sum(nr)
    
  	# A l'altra taula de reports afegeix, per la iteració i cadena avaluada:
  ### En aquest punt aquesta variable correspon als reads que no s'han assignat a la cadena FW!!
  	pr.res[k,1] <- length(seqs)
  	# Reads que s'han pogut associar a la cadena FW
    pr.res[k,2] <- sum(flags)
    # La variable shorts és 0 osigui que el sumatori també ho és
      # Hauria de correspondre als reads que no cobreixen l'amplicó sencer
    pr.res[k,3] <- sum(shorts) #### 
    # Sumatori del total de reads en el MID avaluat (mostra) que s'han assignat
    pr.res[k,4] <- sum(nr)
  
 ### Coincidències cadena reverse
  # Aquestes seqüències corresponen a la cadena reverse, és a dir, presenten el primer RV en la regió 5' i la reversa
  # complementària del primer FW a la regió 3'
    # Cal tenir en compte que la variable seqs ara està actualitzada amb les no assignades a cadena FW!
    # Guarda la seq del primer RV específic de la regió avaluada
    pr.dn <- primers$Primer.RV[ ipr ]
    # Busca la seq del primer RV en la regió 5' (posicions 1-100 definides al fitxer de paràmetres)
    # de les seqüències del fitxer .fna del MID avaluat. max.prdif defineix el màxim de mismatches permesos
    # Guarda totes les coincidències amb les posicions inicial, final, i longitud
    dn.matches <- vmatchPattern(pattern=pr.dn,
                            subject=subseq(seqs,start=target.io,end=target.in),
                            max.mismatch=max.prdif,fixed=FALSE)
    # Resta 1 a la posició inicial de cerca del primer (per defecte, 1)
    delta <- target.io-1
    # Aplica la funció que està definida al fitxer principal
    # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer RV 
    flags <- elementLengths(dn.matches)>=1
    # Torna a buidar aquestes variables per guardar les noves dades
    nr <- integer()
    shorts <- integer()
    # Genera un fitxer .fna, el nom del qual està format per l'ID del pacient, la regió (amplicó) avaluada,
    # i la consecució PrRV.fna (primer reverse)
    dn.flnm <- paste(samples$Patient.ID[jj],primers$Ampl.Nm[ipr],"PrRV.fna",
	                 sep=".")
    
	  seqs.dn <- "" # Coincidències cadena RV
    if(sum(flags)) # Sumatori de seqüències amb 1 o més coincidències amb la seq del primer
    { # 'startIndex()' retorna una llista amb les posicions inicials de les coincidències del patró cercat
      # Cal tenir cura, ja que si no hi ha coincidència la llista guarda un valor NULL
      # Aplicar el condicional flags sobre la llista de posicions inicials permet eliminar els valors NULL
      # Guarda en una variable tots els valors de posició inicial del primer RV sobre les seqüències
      pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta
      
      ### Guarda les seqs amb coincidències del primer RV      
      seqs.dn <- seqs[flags]
      ### Actualitza el nº de seqs amb les que no presentaven coincidències
      seqs <- seqs[!flags]
      
      ### Retalla el primer RV de la regió 5' 
      # Suma a la posició inicial del primer RV en la seq la seva longitud
      st <- pos + (primers$RV.pos[ipr]-primers$RV.tpos[ipr])
      # Retalla les seqüències des de la posició on acaba el primer RV fins al final 
      seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))  
  	  # Guarda només les seqüències amb longitud major a 105 (la posició final on busquem el primer +5)
      seqs.dn <- seqs.dn[width(seqs.dn)>target.in+5]
	  
      ### Actualitza les seqüències de cadena RV fent la seva reversa complementària
        # Important perquè ara les tenim en sentit 3'-5' per poder associar-les a les seves cadenes up (forward)!
      seqs.dn <- reverseComplement(seqs.dn)
      
      ### Busca ara el primer FW original a la regió inicial (5'). En lloc de fer la reversa complementària del primer,
        # ho fa de totes les seqüències on s'ha trobat el primer RV. Per tant, cal buscar la seq original del primer FW
      # Resta 5 unitats a la posició inicial de cerca del primer (per defecte 1) i calcula el màxim entre aquest valor i 1
      io <- max(target.io-5,1)
      # Calcula la diferència entre aquest màxim i la posició inicial de la seqüència
      delta <- io-1
      # Busca la seq del primer FW en la regió 5' (100 primeres posicions) de les seqüències on
      # ja s'ha retallat el primer RV. max.prdif defineix el màxim de mismatches permesos
      # Guarda totes les coincidències amb les posicions inicial, final, i longitud
      dn.matches <- vmatchPattern(pattern=pr.up,
                        subject=subseq(seqs.dn,start=io,end=target.in+5),
                        max.mismatch=max.prdif,fixed=FALSE)
      
      # Indica quines seqüències han trobat 1 o més coincidències amb la seq del primer RV
      flags <- elementLengths(dn.matches)>=1

      ### Els reads que passen aquí tenen l'amplicó sencer, ja que han coincidit els dos primers
      if(sum(flags))
	    { ### Retalla el primer FW (ara a 5') per deixar l'amplicó net
        seqs.dn <- seqs.dn[flags]
        # Torna a fer el procés d'abans i guarda tots els valors de posició inicial del primer FW 
        # sobre les seqüències
        pos <- sapply(startIndex(dn.matches)[flags],function(x) x[[1]])+delta
        # Suma a la posició inicial del primer FW en la seq la seva longitud
        st <- pos + (primers$FW.tpos[ipr]-primers$FW.pos[ipr])
        # Retalla les seqüències des de la posició on acaba el primer FW fins al final 
        seqs.dn <- subseq(seqs.dn,start=st,end=width(seqs.dn))
		
        ### Colapsa les seqüències 'down' en haplotips + freqüències
        # Ordena en ordre descendent les freqüències (calculades amb 'table()') de cada haplotip
          # Les seqüències retallades es transformen en caràcters, i a l'aplicar la funció 'table()', es comparen
          # entre elles, per tant totes les seqs idèntiques s'agrupen juntes. El nº de seqs idèntiques defineix la freq
        sqtbl <- sort(table(as.character(seqs.dn)),decreasing=TRUE)
        # Guarda les seqüències en una variable
        bseqs <- names(sqtbl)                
        # Els noms d'aquesta variable seran nombres de l'1 al total de seqs disponibles
        names(bseqs) <- 1:length(bseqs)
        # Guarda les freqüències dels haplotips (com a nombres enters) ordenades de les seqs
        nr <- as.integer(sqtbl)
        # Aplica la funció del principi sobre les seqüències, les seves freqüències i el fitxer .fna del 
        # pacient al que corresponen aquestes seqs (on es guardaran els resultats)
        lst <- SaveAllHaplotypes(bseqs,nr,dn.flnm)
        # Guarda les freqüències dels haplotips obtinguts
        nr <- lst$nr

        cat("\nReverse seqs, table of read lengths (over 10 rd)\n")
        # 'tapply()' en aquest cas calcula el sumatori (sum) de les freqüències (nº reads) segons les diferents
        # longituds de seqüència
        tbl.len <- tapply(nr,nchar(bseqs),sum)
        # Filtra les longituds amb més de 10 reads (per això el títol amb la funció 'cat()') 
        print(tbl.len[tbl.len>=10])
        # Gràfic histograma en forma de línies verticals representant les diferents longituds en X i les freqüències 
        # (nº reads) en Y, per les seqs reverse
        plot(as.integer(names(tbl.len)),tbl.len,type="h",
             xlab="Read length",ylab="# reads")
        title(main=paste(tt," Str RV"))			 
	  }
	}  
    # Actualitza de nou el valor k per a la propera iteració (mostra) que tornarà a començar per les cadenes up
    k <- k+1
    # Guarda tots els resultats igual que en el cas de les cadenes up (forward) excepte columna $Len
  	FlTbl$File.Name[k] <- dn.flnm
  	FlTbl$Pat.ID[k] <- samples$Patient.ID[jj]
  	FlTbl$Ampl.Nm[k] <- primers$Ampl.Nm[ipr]
  	FlTbl$Pr.ID[k] <- ipr
  	FlTbl$Str[k] <- "rv"
  	# Després de fer la reversa complementària de les cadenes dn (reverse), ambdues cadenes
  	# es troben el mateix sentit, per tant la posició inicial és la mateixa
  	FlTbl$Pos[k] <- primers$FW.tpos[ipr] 
  	# La variable és 0 i no s'actualitza el seu valor
  	  # Si es volgués canviar igual que per les cadenes forward, es podria afegir:
  	  # mean(width(seqs.dn))
  	FlTbl$Len[k] <- trim.len #### 
    FlTbl$Reads[k] <- sum(nr)				   
    FlTbl$Hpls[k] <- length(nr)				   
    
    # Afegeix el nº de reads (del MID concret) idenficats amb la cadena FW al total del pool al qual correspon	  
  	p.ok[i] <- p.ok[i] + sum(nr)
  	# Guarda tots els resultats igual que en el cas de les cadenes up (forward) a l'altra taula (nº de reads per pas)
  	pr.res[k,1] <- length(seqs) #### 
    pr.res[k,2] <- sum(flags)
    pr.res[k,3] <- sum(shorts) #### 
    pr.res[k,4] <- sum(nr)
    
  } # Fi del loop sobre totes les mostres corresponents al pool avaluat
} # Fi del loop sobre els dos pools avaluats
sink()
dev.off()

### Sincronització d'estructures
# Guarda la columna de la taula de resultats amb els identificadors dels primers majors a 0 
# (en aquest cas els IDs son 1 o 2, per tant es guarden tots)
fl <- FlTbl$Pr.ID>0
# Guarda totes les entrades de la taula que compleixin la condició anterior (en aquest cas totes)
FlTbl <- FlTbl[fl,]
# Guarda també les entrades de la taula amb els nº de reads que compleixen la primera condició
pr.res <- pr.res[fl,]
# Guarda tots els identificadors dels pacients concatenas amb la regió del HBV avaluada
anms <- paste(FlTbl$Pat.ID,FlTbl$Ampl.Nm,sep=".")
# Assigna a la taula de mostres (variable samples) els noms de les files, que corresponen de nou als
# identificadors dels pacients amb la regió avaluada (tot i que en aquest cas extrau les dades de la taula de mostres)
rownames(samples) <- paste(samples$Patient.ID,samples$Primer.ID,sep=".")

### Gràfics de resultats
library(RColorBrewer)
# Generació de les paletes de colors
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")

# Guarda en dues variables diferents les entrades de la taula de resultats corresponents a les cadenes forward i reverse
fw.idx <- which(FlTbl$Str=="fw")
rv.idx <- which(FlTbl$Str=="rv")
# Genera un data frame amb dades de les dues taules de resultats que inclou, només per les cadenes forward:
# ID dels pacients, regió amplificada, total de reads de la mostra, dades de la variable shorts (no s'ha definit, és 0),
# reads associats a cadascuna de les dues cadenes (per separat) i sumatori del total de reads en el MID avaluat (mostra) que s'han assignat
mprres <- data.frame(PatID=FlTbl$Pat.ID[fw.idx],
                     PrimerID=FlTbl$Ampl.Nm[fw.idx],
                     # Revisar variable Treads: correspon als reads que no han coincidit per la cadena FW
                     Treads=pr.res[fw.idx,1], #### Perquè fos el total de reads caldria calcular length(seqs) 
                                                 # just després de carregar el fitxer de la carpeta splits
					 Shorts=pr.res[fw.idx,3]+pr.res[rv.idx,3], # En aquest cas és 0 perquè la variable shorts no s'actualitza
                     FW.match=pr.res[fw.idx,4],
                     RV.match=pr.res[rv.idx,4],
					 Fn.reads=pr.res[fw.idx,4]+pr.res[rv.idx,4],
                     stringsAsFactors=FALSE)
# Canvia el nom de la variable
mres <- mprres
# Calcula el sumatori dels reads assignats a cadena forward o reverse en funció dels pacients, és a dir,
# suma els reads assignats per a les dues regions del HBV avaluades
# 'tapply()' calcula el sumatori (sum) del nº reads segons els pacients
T.reads <- apply(mres[,5:6],2, function(x)
              tapply(x,mres$PatID,sum))

# Aquest condicional només s'aplica si a la taula només hi ha un sol pacient
if(length(unique(mres$PatID))==1)
{ # En aquest cas, es genera una matriu d'una sola fila amb els valors de reads assignats a les dues cadenes,
  # indicant el pacient i les columnes de la taula T.reads
  x <- matrix(T.reads,nrow=1)
  rownames(x) <- mres$PatID[1]
  colnames(x) <- names(T.reads)
  T.reads <- x
}  

# Genera el fitxer .pdf que es guardarà a la carpeta reports, indicant el nom del projecte definit al fitxer global
pdf.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.pdf",sep="_")
pdf(file.path(repDir,pdf.flnm),paper="a4",width=6,height=11)
par(mfrow=c(2,1),mar=c(7,4,4,2)+0.1)

# Defineix el límit de l'eix Y del gràfic a partir del sumatori de reads assignats per pacient (no segons el pool)
ymx <- max(rowSums(T.reads))*1.2
# Gràfic de barres amb la trasposada de la taula T.reads, per representar els reads assignats a la cadena forward 
# i reverse per cada pacient, independentment de la regió del HBV avaluada
# En aquest gràfic es representa una única barra per pacient i d'indica amb colors diferents les assignacions up i down
barplot(t(T.reads),col=pal2[1:2],las=2,ylim=c(0,ymx))
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
title(main="Primer matches by patient (# reads)")

# Genera un altre gràfic amb les mateixes dades, però aquest cop es representen dues barres per pacient, una per
# cada cadena forward o reverse, i diferenciant les dues regions del HBV avaluades
res.mat <- mres[,5:6]
ymx <- max(res.mat)*1.2
# També es defineixen els noms de l'eix X com el ID dels pacients i la regió HBV avaluada
nms <- paste(mres$PatID,mres$PrimerID)
bp <- barplot(t(res.mat),beside=TRUE,border=pal1[1:2],col=pal2[1:2],
              ylim=c(0,ymx),xaxt="n")
axis(side=1,at=apply(bp,2,mean),nms,cex.axis=0.6,las=2)
abline(h=0)		  
title(main="Primer matches (# reads)")
legend("top",horiz=TRUE,fill=pal2[1:2],legend=c("up","dn"))
dev.off()

# Calcula el rendiment d'aquest pas, dividint els reads que s'han pogut assignar a alguna de les cadenes 
# entre el total de reads de cada pool o regió (sumatori de tots els MIDS de cada pool)
yield <- p.ok/p.cv*100
# Guarda un data frame amb 3 columnes indicant el total de reads assignats a algun MID, el total de reads assignats
# als primers (cadenes forward o reverse) i el rendiment calculat, tot això per cada regió avaluada
PoolTbl <- data.frame(MIDReads=p.cv,PrimerReads=p.ok,Pct=yield)

# Genera el fitxer .txt que es guardarà a la carpeta reports, indicant el nom del projecte
txt.flnm <- paste(proj.nm,"SplitByPrimersOnFlash.txt",sep="_")
sink(file.path(repDir,txt.flnm))
# Afegeix la taula on es registren, per cada MID avaluat, els reads assignats a forward i reverse 
cat("\nTable of reads identified by primer\n\n")
print(mprres)
cat("\n")
# També afegeix la taula on s'indica la longitud dels reads i els haplotips detectats per cadena (treient el nom del fitxer)
print(FlTbl[,-1])
# Inclou la taula on es mostren els reads assignats per pacient, sumant les dues cadenes
cat("\nTotal reads identified by patient\n\n")
print(T.reads)
# Mostra el rendiment d'aquest pas: els reads totals per pool (tots els MIDS de cada pool) i els que s'han assignat
cat("\nYield by pool\n\n")
print(PoolTbl)
sink()

### Guarda més taules de resultats
# FlTbl <- FlTbl[FlTbl$Reads>0, ]  # Erase rows of null files -> Aquest codi no s'executa
# Guarda les dades en un fitxer .RData
save(FlTbl,PoolTbl,file=file.path(repDir,"SplittedReadsFileTable.RData"))
# Guarda un altre fitxer .txt amb les dades incloses a cada fitxer .fna generat a la carpeta trim
sink(file.path(repDir,"SplittedReadsFileTable.txt"))
print(FlTbl)
sink()


### Genera els gràfics de barres anteriors però en un full A4 horitzontal
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

# També genera dos gràfics boxplot per representar els reads assignats a la cadena forward i reverse per pool
par(mfrow=c(1,2))
# Defineix el límit màxim de l'eix Y
ymx <- max(c(res.mat[,1],res.mat[,2]))
# Defineix els noms de les files de la taula de primers, amb la regió avaluada
rownames(primers) <- primers$Ampl.Nm
# Guarda el nom de la regió avaluada de totes les entrades classificades a la cadena forward
reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="fw"],"Region"]
# Genera el boxplot dels reads assignats a la cadena forward en funció de la regió o pool
boxplot(res.mat[,1]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,1],
       pch="+",cex=0.8)
title(main="Primers identified on forward reads")	   

# Guarda el nom de la regió avaluada de totes les entrades classificades a la cadena reverse
reg <- primers[FlTbl$Ampl.Nm[FlTbl$Str=="rv"],"Region"]
# Genera el boxplot dels reads assignats a la cadena reverse en funció de la regió o pool
boxplot(res.mat[,2]~reg,border="gray",outline=FALSE,
        xlab="",ylab="# of reads",ylim=c(0,ymx))
points(jitter(as.integer(factor(reg)),a=0.15),res.mat[,2],
       pch="+",cex=0.8)
title(main="Primers identified on reverse reads")	   

dev.off()
