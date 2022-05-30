
################################################
###  DISCRIMINANT ANALYSIS BASED ON DISTANCES
###
###  See: Cuadras2008-Cat.pdf pgs 194-195
################################################

### gevoar: estimate of geometric var
###
### Parameters:
###
###   D: distance matrix
###
### OBS: dim(D) = n, where n is the number of individuals of a population
###
geovar <- function(D)
{  if(nrow(D)<2) return(0)
   n <- dim(D)[1]
   D2 <- D ^ 2
   V <- (sum(D2)) / (2 * (n^2))
   return(V)
}

### phi: computes first addend in the estimate of proximity function
###
### Parameters:
###
###   d: vector of distances from one individual to others, that is 
###           d(w, x_1), d(w, x_2),..., d(w, x_n)
###
phi <- function(d)
{  n <- length(d)
   d2 <- d ^ 2
   addend <- sum(d2) / n
   return(addend)
}

### DBdist: computes the nearest cluster to the given sample 
###
### Parameters:
###
###       grpDist: Distance among reference sequences
###            hr: type or subtype for each reference sequence 
###         oDist: distance individual to be classified to the reference sequences.
###      my.names: cluster names, defaults to its index.
###
### Return value: a vector with the distances to each cluster, the index
###       of the nearest cluster and its name.
###
DBrule <- function(grpDist, hr, oDist, my.names=NULL)
{
   g <- max(hr)
   grpDist <- as.matrix(grpDist)
   PHI2 <- vector(mode="numeric",length=g)
   names(PHI2) <- paste("Phi2",1:g,sep=".")
   if(!is.null(my.names))    names(PHI2) <- paste("Phi2",my.names,sep=".")
   for(i in 1:g)
   { D <- grpDist[hr==i,hr==i,drop=FALSE]
     V <- geovar(D)
     PHI2[i] <- phi(oDist[hr==i]) - V
   }
   idx <-  match(min(PHI2), PHI2)
   clustName <- ifelse(!is.null(my.names),my.names[idx],idx)
   output <- list(Phi2=PHI2, DB.rule=idx, Type=clustName)
   return(output)   
}

