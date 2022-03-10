

load(file.path(resultsDir,"5pHBx_Classif.RData"))
x.gtbl <- gtbl
x.gtbl2 <- gtbl2
x.df.gtbl <- df.gtbl
x.df.gtbl2 <- df.gtbl2

load(file.path(resultsDir,"PreS1_Classif.RData"))
s.gtbl <- gtbl
s.gtbl2 <- gtbl2
s.df.gtbl <- df.gtbl
s.df.gtbl2 <- df.gtbl2

format.ln <- function(v)
{ v <- v[v>0]
  ln <- paste(names(v)," ",v,"%,",sep="")
  ln <- substr(ln,1,nchar(ln)-1)
  if(length(ln)>1)
    ln <- paste(ln,collapse="; ")
  return(ln)
}  

format.ln.rd <- function(v)
{ v <- v[v>0]
  ln <- paste(names(v)," ",v,",",sep="")
  ln <- substr(ln,1,nchar(ln)-1)
  if(length(ln)>1)
    ln <- paste(ln,collapse="; ")
  return(ln)
}  

IDs <- sort(union(rownames(x.gtbl),rownames(s.gtbl)))
n <- length(IDs)
pct <- data.frame(x5pX=rep("",n),x.flag=rep("",n),
                  PreS1=rep("",n),s.flag=rep("",n),
                  stringsAsFactors=FALSE)
rds <- data.frame(x5pX=rep("",n),x.flag=rep("",n),
                  PreS1=rep("",n),s.flag=rep("",n),
                  stringsAsFactors=FALSE)
rownames(pct) <- rownames(rds) <- IDs				   
for(id in IDs)
{ if(id %in% rownames(x.gtbl))
  { pct[id,"x5pX"] <- format.ln(x.gtbl2[id,])
    rds[id,"x5pX"] <- format.ln.rd(x.gtbl[id,])
	pct[id,"x.flag"] <- x.df.gtbl2[id,"flag"]
    rds[id,"x.flag"] <- x.df.gtbl[id,"flag"]
  }
  if(id %in% rownames(s.gtbl))
  { pct[id,"PreS1"] <- format.ln(s.gtbl2[id,]) 
    rds[id,"PreS1"] <- format.ln.rd(s.gtbl[id,]) 
	pct[id,"s.flag"] <- s.df.gtbl2[id,"flag"]
    rds[id,"s.flag"] <- s.df.gtbl[id,"flag"]
  }
}

sink(file.path(resultsDir,"HBV.Genotyping.Rprt.txt"))
cat("\nGenotyping results by 5\'X and PreS1 amplicons\n")
cat("\nAs percentage:\n")
print(pct)
cat("\nIn reads number:\n")
print(rds)
sink()


