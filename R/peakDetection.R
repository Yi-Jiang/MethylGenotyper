
#' Estimate mode location for each probe
#' 
#' @param x Matrix of beta or RAI values. Row names must be supplied.
#' @param bw band width.
#' @param minDens Minimum density for a valid peak.
#' @param maxProp_antimode If the proportion of antimode density on mode density exceed this threshold, the two peaks will be merged.
#' @param cpu Number of CPU.
#' @return A data frame of mode locations.
#' @export
getMod <- function(x, bw=0.04, minDens=0.01, maxProp_antimode=0.5, cpu=1){
  cl<- makeCluster(cpu)
  registerDoParallel(cl) 
  modRes <- foreach(cpg=rownames(x), .packages=c("tidyverse","multimode"), 
                    .export=c("findCentralFromTwoPeaks", "findCentralFromTwoPeaks")) %dopar% {
    tryCatch({
      # Get mode location and density
      nMod <- nmodes(x[cpg,], bw, lowsup=0, uppsup=1)
      loc <- locmodes(x[cpg,], mod0=nMod)
      if(length(loc$locations) < 2){
        print(paste0("Escape ", cpg, " as <2 valid peaks detected."))
        return(c(CpG=cpg, nmod=length(loc$locations), loc_pass=FALSE, loc0=NA, loc1=NA, loc2=NA))
      }
      idx <- seq(1, length(loc$locations), 2)
      idx2 <- seq(2, length(loc$locations), 2)
      modes <- data.frame(loc=loc$locations[idx], dens=loc$fvalue[idx])
      antimodes <- data.frame(loc=loc$locations[idx2], dens=loc$fvalue[idx2])
      
      # Remove modes with low density. For the two nearby antimodes, remove the higher one.
      lowDens <- modes$dens <= minDens
      modes <- modes[!lowDens,]
      antimodes <- antimodes[!lowDens[-1],]
      if(lowDens[1]){antimodes <- antimodes[-1,]}
      
      # # Merge peaks if antimode density > 1/2 of mode density
      # # Not dealing well with low density peaks, like rs556 (cg14655569).
      # antimodes$pLeft <- antimodes$dens / modes$dens[-nrow(modes)]
      # antimodes$pRight <- antimodes$dens / modes$dens[-1]
      # if(max(antimodes[,c("pLeft", "pRight")]) > maxProp_antimode){
      #   for(i in nrow(antimodes):1){
      #     if(max(antimodes[i, c("pLeft", "pRight")]) > maxProp_antimode){
      #       antimodes <- antimodes[-i,]
      #       if(modes[i, "dens"] > modes[i+1, "dens"]){
      #         modes <- modes[-(i+1),]
      #       }else{
      #         modes <- modes[-i,]
      #       }
      #     }
      #   }
      # }

      # Detect the central mode
      if(nrow(modes)==2){
        loc012 <- findCentralFromTwoPeaks(modes$loc)
      }else if(nrow(modes)==3){
        loc012 <- sort(modes$loc)
      }else if(nrow(modes)>3){ # remove the lowest peaks as they are mostly noise
        while(nrow(modes)>3){
          lowest <- which.min(modes$dens)
          if(lowest==1){
            antimodes <- antimodes[-1,]
          }else{
            antimodes <- antimodes[-(lowest - 1),]
          }
          modes <- modes[-lowest,]
        }
        loc012 <- sort(modes$loc)
      }else{
        print(paste0("Escape ", cpg, " as <2 valid peaks detected."))
        return(c(CpG=cpg, nmod=nrow(modes), loc_pass=FALSE, loc0=NA, loc1=NA, loc2=NA))
      }
      #if(loc012[2]>0.3 & loc012[2]<0.7){
        #if(all(c(loc012[1]<.3, loc012[3]>.7), na.rm=TRUE)){
          loc_pass=TRUE
        #}else{loc_pass=FALSE}
      #}else{loc_pass=FALSE}
      return(c(CpG=cpg, nmod=nrow(modes), loc_pass=loc_pass, loc0=loc012[1], loc1=loc012[2], loc2=loc012[3]))
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


