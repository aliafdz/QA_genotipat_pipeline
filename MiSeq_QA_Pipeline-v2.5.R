
#####################################################
###     MiSeq - Quality assessment pipeline       ###
###             and Q30 filter                    ###
#####################################################

library(Biostrings)
library(ShortRead)

## Definim els directoris (carpetes) indicant la ruta:

# Directori on trobem els arxius de codi que s'aniran cridant. 
codeDir <- "./R" 
# Carpeta data on trobem els primers, pacients, mids, seq de referència, codi genètic, etc. 
dataDir <- "./data"
# On es depositen els arxius fastq.gz
runDir <- "./run"
# Crec que aquests dos només incorporen alguns arxius que generem en els passos del pipeline. 
flashDir <- "./flash"
flashFiltDir <- "./flashFilt"
# Arxius resultants del pipeline a nivell de qualitat. 
repDir <- "./reports"

min.len <- 200   # Longitud mínima per considerar una seqüència
min.ov <- 20     # Mínim solapament (en nt) entre R1 and R2
max.ov <- 300    # Màxim solapament (en nt) entre R1 and R2
err.lv <- 0.10   # Fracció de diferències acceptades en el solapament
flash <- "C:/FLASH/flash.exe" # Carpeta on es troba l'executable FLASH.

# Indicador de les variables min, max i error level
flash.opts <- paste("-m",min.ov,"-M",max.ov,"-x",err.lv)  

### Chunck size to be used by FastqStreamer() --> Buscar info
# Nombre de registres successius a retornar a cada rendiment (yield)
chunck.sz <- 1.e6

tm <- integer(10) # Vector amb 10 nombres 0

###  Control d'objectes necessaris
par.nms <- c(ls(),"par.nms")

###  Per controlar l'espai ocupat per cada objecte en memòria:
###  sort(sapply(ls(),function(onm) object.size(eval(get(onm)))),
###       decreasing=TRUE)
###  Veure també gc() i memory.profile()

###  Alliberar memòria innecessària després de cada pas
###      rm(list=setdiff(ls(),par.nms))   
###  notar que dins d'una funció no va s'espera.

# cat("\nRunning quality assessment on R1 and R2 fastq files\n")
# print( (tm[1] = system.time( source("./R/R1R2_LoqQ2N_pl-v1.27.R") )[3]) ) 

cat("\nRunning FLASH to extend reads\n")
print( (tm[2] = system.time( source("./R/R1R2_to_FLASH_pl-v2.R") )[3]) ) 
rm(list=setdiff(ls(),par.nms))

cat("\nRunning QC by position, R1, R2 and FLASH fastq files\n")
print( (tm[3] = system.time( source("./R/PoolQCbyPos-v2.2.R") )[3]) ) 
rm(list=setdiff(ls(),par.nms))

cat("\nRunning QC by read on FLASH fastq files\n")
print( (tm[4] = system.time( source("./R/PoolQCbyRead-v1.15.R") )[3]) ) 
rm(list=setdiff(ls(),par.nms))

cat("\nRunning length peaks on FLASH fastq files\n")
print( (tm[5] = system.time( source("./R/LenPeaksByPool-v1.R") )[3]) ) 
rm(list=setdiff(ls(),par.nms))

###  Accept max of 5% bases below Q30 by read
ThrQ30 <- 0.05
cat("\nFiltering haplotypes by Q30\n")
print( (tm[6] = system.time( source("./R/FiltByQ30-v1.1.R") )[3]) ) 

cat("\nRunning QC by position after filtering fastq files\n")
print( (tm[7] = system.time( source("./R/PoolFiltQCbyPos-v2.2.R") )[3]) ) 

###  End of pipeline
cat("\nEnd of MiSeq quality assessment pipeline")
cat("\nElapsed globally: ",sum(tm),"\" (",round(sum(tm)/60,2),"\', ",
    round(sum(tm)/3600,2),"h)\n",sep="")
	
#---------------------------------------------------------------------------#

 # tm
 # [1]    0.00 693.75 849.51 287.80  73.58 396.46  78.29
 
 # End of MiSeq quality assessment pipeline
 #   Elapsed globally: 2379.39" (39.66', 0.66h)
