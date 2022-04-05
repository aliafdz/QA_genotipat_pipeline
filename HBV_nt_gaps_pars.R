
############################################################
###  PARAMETRES DEL PIPE-LINE DE InDels i Genotypat de HBV
############################################################

stopifnot(require(Biostrings))
stopifnot(require(stringr))
stopifnot(require(RColorBrewer))

dataDir <- "./data" # Carpeta amb arxius de mostres, primers i MIDS
codeDir <- "./R"    # Arxius de codi R
runDir <- "./run"   # Arxius de seqüenciació fastq.gz
flashDir <- "./flashFilt" # Arxius filtrats per Q30. Important!! Es redefineix el directori respecte QA 
splitDir <- "./splits" # Seqüències demultiplexades per MID
trimDir <- "./trim"    # Seqüències demultiplexades per primers i retallades
filtDir <- "./filt" # Aquest directori no es fa servir!
joinDir <- "./join" # Aquest directori no es fa servir!
repDir <- "./reports" # Carpeta on es guarden els resultats del pipeline
resultsDir <- resDir <- "./results" # Guarda els resultats més rellevants
exportDir <- expDir <- "./export" # Aquest directori no es fa servir!
tempDir <- "./tmp" # Arxius temporals generats amb MUSCLE
ntDir <- "./nt" # Aquest directori no es fa servir!
aaDir <- "./aa" # Aquest directori no es fa servir!

data.Dir <- repDir 
desc.Dir <- dataDir
mach.Dir <- "./MACH" # Seqüències FW i RV aliniades i intersecció d'haplotips
tmp.Dir <- "./tmp"

#-----------------------------------------------------------------------#

###  Posició esperada del MID (general, adaptador o no)
mid.start <- 1
mid.end <- 40

#-----------------------------------------------------------------------#

##  TRIM PRIMERS
## Definim paràmetres per separació de reads:
pmm.mx <- 3        ## Nombre màxim de mismatch en el primer específic
max.prdif <- pmm.mx
min.len <- 200     ## Longitud mínima per considerar una seqüència
###  Depèn que hi hagi adaptadors, MIDs i/o M13
target.io <- 1     #  1   10   25   50
target.in <- 100   # 30   55   55   80

#-----------------------------------------------------------------------#

## Definim paràmetres per la intersecció d'haplotips

## Mínima longitud per entrar en la intersecció
min.seq.len <- 150
###  Si cal filtrat per mínim nombre de reads en la lectura del fasta
min.rd <-   1      
###  Min d'abundància per entrar en l'aliniament múltiple (%)
a.cut  <- 0.2

#-----------------------------------------------------------------------#
