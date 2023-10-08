
#' Estimate mode location for each probe
#' 
#' @param x Matrix of beta or RAI values. Row names must be supplied.
#' @param minDens Minimum density for a valid peak.
#' @param cpu Number of CPU.
#' @return A data frame of mode locations.
#' @export
getMod <- function(x, minDens=0.01, cpu=1){
  print(paste(Sys.time(), "Running mode test."))
  cl<- makeCluster(cpu)
  registerDoParallel(cl) 
  modRes <- foreach(cpg=rownames(x), .packages=c("tidyverse","multimode"), 
                    .export=c("findCentralFromTwoPeaks", "findCentralFromTwoPeaks")) %dopar% {
    tryCatch({
      # Estimate number of modes for each probe
      nMod <- nmodes(x[cpg,], 0.04, lowsup=0, uppsup=1)
      loc <- locmodes(x[cpg,], mod0=nMod)
      idx <- seq(1, length(loc$locations), 2)
      modes <- data.frame(loc=loc$locations[idx], dens=loc$fvalue[idx])
      modes <- modes[modes$dens>minDens,]
      if(nrow(modes)==2){
        loc012 <- findCentralFromTwoPeaks(modes$loc)
      }else if(nrow(modes)==3){
        loc012 <- sort(modes$loc)
      }else if(nrow(modes)>3){
        while(nrow(modes)>3){
          modes <- rmSmallerFromClosestPeaks(modes)
        }
        loc012 <- sort(modes$loc)
      }else{
        print(paste0("Escape ", cpg, " as <2 valid peaks detected."))
        return(c(CpG=cpg, nmod=nrow(modes), loc_pass=FALSE, loc0=NA, loc1=NA, loc2=NA))
      }
      if(loc012[2]>0.3 & loc012[2]<0.7){loc_pass=TRUE}else{loc_pass=FALSE}
      return(c(CpG=cpg, nmod=nrow(modes), loc_pass=loc_pass, loc012[1], loc012[2], loc012[3]))
    }, error = function(e) return(paste0("Escape ", cpg, " with error: ", e)))
  }
  stopImplicitCluster()
  stopCluster(cl)
  modRes <- as.data.frame(do.call(rbind, modRes))
  modRes$nmod <- as.numeric(modRes$nmod)
  modRes$loc0 <- as.numeric(modRes$loc0)
  modRes$loc1 <- as.numeric(modRes$loc1)
  modRes$loc2 <- as.numeric(modRes$loc2)
  rownames(modRes) <- modRes$CpG
  modRes
}

#' Find the central peak from a distribution of two peaks
#' 
#' @param locations A vector of two, indicating the locations of the two peaks.
#' @return A vector of three, with the central one being the location of the central peak.
#' @export
findCentralFromTwoPeaks <- function(locations){
  if(length(locations)!=2){
    print("Error: findCentralFromTwoPeaks() only accept a vector of two as input.")
    return(NA)
  }
  locations <- c(NA, locations, NA)
  distance <- abs(locations - 0.5)
  if(distance[2] < distance[3]){
    loc012 <- locations[1:3]
  }else{
    loc012 <- locations[2:4]
  }
  names(loc012) <- paste0("loc", 0:2)
  loc012
}

#' Find the closest peaks and remove the smaller one
#' 
#' @param modes A matrix of two columns (col names: loc, dens) and at least two rows. 
#' @return The same matrix with one row deleted.
#' @export
rmSmallerFromClosestPeaks <- function(modes){
  if(intersect(c("loc", "dens"), colnames(modes))!=2 | nrow(modes)<2){
    print("Error: rmSmallerFromClosestPeaks() requres a matrix of two columns and at least two rows as input.")
    return(NA)
  }
  closest <- which.min(diff(modes$loc))
  if(which.min(modes[c(closest, closest+1), "dens"])==1){
    modes <- modes[-closest, ]
  }else{
    modes <- modes[-(closest+1), ]
  }
  modes
}


