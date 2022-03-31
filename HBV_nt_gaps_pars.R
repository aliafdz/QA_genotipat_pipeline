
############################################################
###  PARAMETRES DEL PIPE-LINE DE InDels i Genotypat de HBV
############################################################

stopifnot(require(Biostrings))
stopifnot(require(stringr))
stopifnot(require(RColorBrewer))

dataDir <- "./data"
codeDir <- "./R"
runDir <- "./run"
flashDir <- "./flashFilt" # Important!! Es redefineix el directori respecte QA
splitDir <- "./splits"
trimDir <- "./trim"
filtDir <- "./filt" # Aquest directori no es fa servir!
joinDir <- "./join"
repDir <- "./reports"
resultsDir <- resDir <- "./results"
exportDir <- expDir <- "./export"
tempDir <- "./tmp"
ntDir <- "./nt"
aaDir <- "./aa"

data.Dir <- repDir
desc.Dir <- dataDir
mach.Dir <- "./MACH"
tmp.Dir <- "./tmp"

#-----------------------------------------------------------------------#

###  Posició esperada del MID (general, adaptador 454 o no)
mid.start <- 1
mid.end <- 40

#-----------------------------------------------------------------------#

##  TRIM PRIMERS
## Definim parametres per separació de reads:
pmm.mx <- 3        ## Nombre màxim de mismatch en el primer específic
max.prdif <- pmm.mx
min.len <- 200     ## Longitud mínima per considerar una seqüència
###  Depèn que hi hagi adaptadors 454, MIDs i/o M13
target.io <- 1     #  1   10   25   50
target.in <- 100   # 30   55   55   80

#-----------------------------------------------------------------------#

## Definim paràmetres per la intersecció d'haplotips

## Mínima longitud per entrar en la intersecció
min.seq.len <- 150
###  Si cal filtrat per mínim nombre de reads en la lectura del fasta
min.rd <-   1      
###  Min d'abundancia per entrar en l'aliniament múltiple (%)
a.cut  <- 0.2

#-----------------------------------------------------------------------#
