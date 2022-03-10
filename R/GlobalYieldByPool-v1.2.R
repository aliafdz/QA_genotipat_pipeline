
###  Llegir la descripció de mostres
samples <- read.table(file.path(desc.Dir,"samples.csv"), sep="\t",
                      header=T,stringsAsFactors=F)
primers <- read.table(file.path(desc.Dir,"primers.csv"), sep="\t",
                      header=T,stringsAsFactors=F)

###  Carregar taules
load(file=file.path(repDir,"FLASH_table.RData"))
load(file=file.path(repDir,"SplittedReadsFileTable.RData"))
rownames(flash.res) <- str_extract(rownames(flash.res),"^[A-Za-z0-9\\.-]+")
load(file=file.path(repDir,"IntersectedReads.RData"))

###  Rendiment per pas de cada pool
#####################################

p.nms <- rownames(PoolTbl)
flash.res <- flash.res[p.nms,]

rds.raw <- flash.res[,"Extended"]+flash.res[,"NoExtd"]
rds.flash <- flash.res[,"Extended"]
rds.fQ30 <- flash.res[,"FiltQ30"]
rds.MID <- PoolTbl[p.nms,"MIDReads"]
rds.dmult <- PoolTbl[p.nms,"PrimerReads"]

frdf.pool <- sapply(1:nrow(frdf), function(i)
               samples$Pool.Nm[ which(samples$Patient.ID==frdf$Pat.ID[i] &
		                              samples$Primer.ID==frdf$Ampl.Nm[i])[1] ])
								  
rds.filt <- tapply(frdf$all,frdf.pool,sum)[p.nms]
rds.ints <- tapply(frdf$Fn.rd,frdf.pool,sum)[p.nms]
gbly <- cbind(raw=rds.raw,flash=rds.flash,fQ30=rds.fQ30,
              MID=rds.MID,primer=rds.dmult,
              filt=rds.filt,ints=rds.ints)
Tgbl <- colSums(gbly)
pct.gbl <- round(Tgbl[-1]/Tgbl[-length(Tgbl)]*100,2)

pdf.flnm3 <- "GlobalYieldBarplots.pdf"
pdf(file.path(repDir,pdf.flnm3),paper="a4",width=6.5,height=10)
par(mfrow=c(2,1))

pal1 <- brewer.pal(8,"Dark2")
pal2 <- brewer.pal(8,"Pastel2")
m <- ncol(gbly)
ymx <- max(gbly)*1.2
bp <- barplot(t(gbly),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,ymx),xaxt="n",cex.axis=0.8,ylab="# reads")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
abline(h=0)	
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
title(main="Yield on pools by analysis step")

bp <- barplot(t(gbly/gbly[,1]*100),beside=TRUE,col=pal2[1:m],border=pal1[1:m],
              ylim=c(0,117),xaxt="n",cex.axis=0.8,ylab="Percentage")
axis(side=1,at=apply(bp,2,mean),rownames(gbly),las=2,cex.axis=0.8)	
grid(nx=NA,ny=NULL)	
abline(h=0)	
legend("top",horiz=TRUE,fill=pal2[1:m],legend=colnames(gbly),cex=0.8)		
title(main="Yield on pools by analysis step")

par(mar=c(5,8,4,6))

bp <- barplot(pct.gbl,col="lavender",border="navy",ylim=c(0,max(pct.gbl)),
        ylab="yield (%)")
text(bp,50,pct.gbl,col="navy",font=2,cex=0.8)
title(main="Global yield by step")		

logReads <- ifelse(frdf$Fn.rd>0,log10(frdf$Fn.rd),0)
boxplot(logReads,border="gray",ylab="log10(# of reads)",outline=FALSE,
        ylim=range(logReads))
points(jitter(rep(1,nrow(frdf)),a=0.10),logReads,pch="+",cex=0.8)
title(main="Final coverage")

dev.off()
			  
txt.flnm2 <- "GlobalYield-SumRprt.txt"
sink(file=file.path(repDir,txt.flnm2))

cat("\n   Global yield by analysis step")
cat("\n===================================\n")
cat("\nIn number of reads:\n\n")
gbly <- rbind(gbly,TOTAL=colSums(gbly))
print(gbly)	
		
cat("\nIn percentage by step:\n\n")
print(round(gbly[,-1]/gbly[,-ncol(gbly)]*100,2))

cat("\nIn percentage referred to raw reads:\n\n")
print(round(gbly/gbly[,1]*100,2))			
sink()
