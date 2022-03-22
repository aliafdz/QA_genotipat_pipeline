
### Llegeix l'estructura de descripció de mostres
samples <- read.table(file.path(desc.Dir,"samples.csv"), sep="\t",
                      header=T,stringsAsFactors=F)

### Llegeix el fitxer amb els primers
primers <- read.table(file.path(desc.Dir,"primers.csv"), sep="\t",
                      header=T,stringsAsFactors=F)

### Carrega les taules de resultats .RData després de FLASH, després del trimming de primers
# i després de fer la intersecció dels haplotips
load(file=file.path(repDir,"FLASH_table.RData")) # Inclou taula flash.res
load(file=file.path(repDir,"SplittedReadsFileTable.RData")) # Inclou taules FlTbl i PoolTbl
load(file=file.path(repDir,"IntersectedReads.RData")) # Inclou taula frdf
# La funció 'str_extract()' permet treure, dels noms de les files de la taula, el terme
# final _S2 o _S1 (s'entén que correspon al terme "^[A-Za-z0-9\\.-]+")
rownames(flash.res) <- str_extract(rownames(flash.res),"^[A-Za-z0-9\\.-]+")


###  Rendiment per pas de cada pool
# Guarda els noms de les files de la taula PoolTbl (noms dels pools o regions avaluades)
p.nms <- rownames(PoolTbl)
# Guarda la taula flash.res, en concret les files corresponents als pools anteriors
# (en aquest cas, tota la taula)
flash.res <- flash.res[p.nms,]

# Guarda en una variable les columnes de la taula flash.res corresponents al total de reads
# estesos i no estesos per FLASH (la suma dels dos son els reads inicials o raw reads) per pool
# Nota: aquest resultat es podria agafar de la columna Raw i no fent la suma
rds.raw <- flash.res[,"Extended"]+flash.res[,"NoExtd"]
# Guarda els reads estesos per FLASH
rds.flash <- flash.res[,"Extended"]
# Guarda els reads restants després del filtrat per Q30
rds.fQ30 <- flash.res[,"FiltQ30"]
# Guarda els reads restants després del demultiplexat per MIDS (total per cada pool)
rds.MID <- PoolTbl[p.nms,"MIDReads"]
# Guarda els reads restants després del demultiplexat per primers
rds.dmult <- PoolTbl[p.nms,"PrimerReads"]

# Guarda el noms dels pools corresponents a cada mostra (?
# which(samples$Patient.ID==frdf$Pat.ID[i] & samples$Primer.ID==frdf$Ampl.Nm[i])[1]
# indica quin identificador de la taula samples correspon al pacient corresponent de la taula frdf
# i també a la regió avaluada segons el valor i
# Guarda el pool corresponent després d'aplicar el condicional per a totes les files de la taula frdf (mostres) 
frdf.pool <- sapply(1:nrow(frdf), function(i)
               samples$Pool.Nm[ which(samples$Patient.ID==frdf$Pat.ID[i] &
		                              samples$Primer.ID==frdf$Ampl.Nm[i])[1] ])

# Guarda el total de reads després del filtre per longitud dels haplotips  		  
rds.filt <- tapply(frdf$all,frdf.pool,sum)[p.nms]
# Guarda el total de reads després de la intersecció entre haplotips
rds.ints <- tapply(frdf$Fn.rd,frdf.pool,sum)[p.nms]

## Genera una taula on s'incloguin tots els resultats guardats anteriorment en columnes
gbly <- cbind(raw=rds.raw,flash=rds.flash,fQ30=rds.fQ30,
              MID=rds.MID,primer=rds.dmult,
              filt=rds.filt,ints=rds.ints)
# Guarda el sumatori dels resultats dels dos pools avaluats
Tgbl <- colSums(gbly)
# Calcula el rendiment (en %) de tots els passos respecte l'anterior (pools sumats)
pct.gbl <- round(Tgbl[-1]/Tgbl[-length(Tgbl)]*100,2)

## Genera el fitxer pdf on es guardaran els gràfics de resultats a la carpeta reports
pdf.flnm3 <- "GlobalYieldBarplots.pdf"
pdf(file.path(repDir,pdf.flnm3),paper="a4",width=6.5,height=10)
par(mfrow=c(2,1))
# Genera les paletes de colors
pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
# Total de columnes de la taula gbly, és a dir els passos totals de l'anàlisis
m <- ncol(gbly)
# Defineix el límit de l'eix Y a partir del màxim de la taula de resultats
ymx <- max(gbly)*1.2
## Gràfic de barres representant el nº de reads de cada pas de l'anàlisis per pool
bp <- barplot(t(gbly),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,ymx),xaxt="n",cex.axis=0.8,ylab="# reads")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
abline(h=0)	
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
title(main="Yield on pools by analysis step")

## Gràfic de barres representant el rendiment (en %) de cada pas de l'anàlisis per pool
bp <- barplot(t(gbly/gbly[,1]*100),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,117),xaxt="n",cex.axis=0.8,ylab="Percentage")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
grid(nx=NA,ny=NULL)	
abline(h=0)	
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
title(main="Yield on pools by analysis step")

## Gràfic de barres representant el rendiment global de cada pas de l'anàlisis (pools sumats)
par(mar=c(5,8,4,6))
bp <- barplot(pct.gbl,col="lavender",border="navy",ylim=c(0,max(pct.gbl)),
        ylab="yield (%)")
text(bp,50,pct.gbl,col="navy",font=2,cex=0.8)
title(main="Global yield by step")		

## Diagrama de caixa representant la cobertura final en funció dels reads de l'últim pas
# 'ifelse()' permet obtenir un resultat en funció del test del primer argument (TRUE o FALSE)
# En aquest cas, si la columna final reads de la taula frdf és major a 0, retorna el logaritme en base 10,
# i si no retorna com a resultat 0. 
logReads <- ifelse(frdf$Fn.rd>0,log10(frdf$Fn.rd),0)
boxplot(logReads,border="gray",ylab="log10(# of reads)",outline=FALSE,
        ylim=range(logReads))
points(jitter(rep(1,nrow(frdf)),a=0.10),logReads,pch="+",cex=0.8)
title(main="Final coverage")

dev.off()


## Genera el fitxer .txt on es guardaran els resultats finals a la carpeta reports
txt.flnm2 <- "GlobalYield-SumRprt.txt"
sink(file=file.path(repDir,txt.flnm2))

cat("\n   Global yield by analysis step")
cat("\n===================================\n")
cat("\nIn number of reads:\n\n")
# Guarda la taula de nº de reads segons el pas de l'anàlisi per pool, i afegeix el sumatori 
# dels dos pools a la tercera fila 
# Nota: es podria substituir per la variable Tgbl
gbly <- rbind(gbly,TOTAL=colSums(gbly))
print(gbly)	
		
cat("\nIn percentage by step:\n\n")
# Guarda els resultats de % dels passos en funció de l'anterior, per pool
print(round(gbly[,-1]/gbly[,-ncol(gbly)]*100,2))

cat("\nIn percentage referred to raw reads:\n\n")
# Guarda els resultats de % dels passos en funció dels raw reads, per pool
print(round(gbly/gbly[,1]*100,2))			
sink()
