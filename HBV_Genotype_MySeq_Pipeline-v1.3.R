
########################################################
###     HBV - GAPS ANALYSIS + GENOTIPAT PIPELINE     ###
########################################################

# Defineix el nom del projecte
proj.nm <- "BQ68"

source("HBV_nt_gaps_pars.R")  # Executa el codi on es defineixen els paràmetres del processament

#source("NewCode2OldRvers.R")  # Per córrer en versions anteriors a la 3.

# La funció 'elementLengths()' està desactualitzada i ara correspon a la funció 'elementNROWS()', 
# d'aquesta manera s'eviten errors en l'execució
if(!exists("elementLengths"))
  if(exists("elementNROWS"))
    { elementLengths <- function(x)
	    return(elementNROWS(x)) # -> Retorna el vector que conté el nombre
                              # de coincidències per cada patró cercat
	}

### Aquesta part està modificada respecte el fitxer original
### Afegim l'executable muscle necessari per alguns arxius del pipeline, així no es donaran errors per no trobar
  # la ruta de l'arxiu al canviar d'ordinador
muscle <- "C:\\Muscle\\muscle3.8.31_i86win32.exe"
# Afegim un test per comprovar que el fitxer existeix i es pot executar correctament
if(!file.exists(muscle)){
  stop("Indica la ruta correcta de l'executable muscle")
}

tm <- integer(10)

cat("\nSplit by MIDs\n")
print( (tm[1] = system.time( 
         source("./R/BigFilesSplitMIDs-v1.27.R") )[3]) ) 

cat("\nTrimming adaptors and primers\n")
print( (tm[1] = system.time( 
         source("./R/RAV_PrimersSplitFromMIDsSplits_WithGaps-pl-v2.30.R") )[3]) ) 

cat("\nConsensus haplotypes by multiple alignment\n")
print( (tm[2] = system.time( source("./R/RawConsHaplosByMA-v1.6.R") )[3]) ) 

cat("\nGlobal yield by step in pools\n")
print( (tm[3] = system.time( source("./R/GlobalYieldByPool-v1.2.R") )[3]) ) 

cat("\nPlot insertions and deletions\n")
print( (tm[4] = system.time( source("./R/InsDelsByMA-v1.2.R") )[3]) ) 

cat("\nCoverage, haplotypes and max difs\n")
print( (tm[5] = system.time( source("./R/NtFastasSummary-v1.2.R") )[3]) ) ###


### Els arxius a partir d'aquí no s'han inclòs en la revisió realitzada en el TFM ###
cat("\nGenotyping by 5pHBx\n")
print( (tm[6] = system.time( source("./R/ClassifHpl-HBV-5pHBx-v10.R") )[3]) ) 

cat("\nGenotyping by PreS1\n")
print( (tm[7] = system.time( source("./R/ClassifHpl-HBV-PreS1-v10.R") )[3]) ) 

cat("\nIntegrate genotyping by both regions\n")
print( (tm[8] = 
   system.time( source("./R/Integrate_5pX_PreS1_Genotyping-v2.R") )[3]) ) 

cat("\nTest possible contaminations\n")
print( (tm[9] = 
   system.time( source("./R/AllHaplotype.Tree&Equivalencies-v2.1.R") )[3]) ) 

###  End of pipeline
cat("\nEnd of MiSeq nucleotide pipeline")
cat("\nElapsed globally: ",sum(tm),"\" (",round(sum(tm)/60,2),"\', ",
    round(sum(tm)/3600,2),"h)\n",sep="")

#---------------------------------------------------------------------------#
