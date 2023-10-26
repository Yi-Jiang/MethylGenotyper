
#' Estimate mode location for each probe
#' 
#' @param x Matrix of beta or RAI values. Row names must be supplied.
#' @param bw band width.
#' @param minDens Minimum density for a valid peak.
#' @param cpu Number of CPU.
#' @return A data frame of mode locations.
#' @export
getMod <- function(x, bw=0.04, minDens=0.001, cpu=1){
  cl<- makeCluster(cpu)
  registerDoParallel(cl) 
  modRes <- foreach(cpg=rownames(x), .packages=c("tidyverse","multimode"), 
                    .export=c("findCentralFromTwoPeaks")) %dopar% {
    tryCatch({
      # Get mode location and density
      nMod <- nmodes(x[cpg,], bw, lowsup=0, uppsup=1)
      loc <- locmodes(x[cpg,], mod0=nMod)
      if(length(loc$locations) < 2){
        print(paste0("Escape ", cpg, " as <2 valid peaks detected."))
        return(c(CpG=cpg, nmod=length(loc$locations), loc_pass=FALSE, loc0=NA, loc1=NA, loc2=NA, dens0=NA, dens1=NA, dens2=NA))
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
      modes <- modes[order(modes$loc),]
      antimodes <- antimodes[order(antimodes$loc),]

      # Detect the central mode
      if(nrow(modes)==2){
        distance <- abs(modes$loc - 0.5)
        if(distance[1] < distance[2]){
          modes <- rbind(NA, modes)
        }else{
          modes <- rbind(modes, NA)
        }
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
      }else if(nrow(modes)<2){
        print(paste0("Escape ", cpg, " as <2 valid peaks detected."))
        return(c(CpG=cpg, nmod=nrow(modes), loc_pass=FALSE, loc0=NA, loc1=NA, loc2=NA, dens0=NA, dens1=NA, dens2=NA))
      }
      loc012 <- modes$loc
      dens012 <- modes$dens
      
      # Position of the central mode in a given range?
      if(loc012[2]>0.3 & loc012[2]<0.7){
          loc_pass=TRUE
      }else{
        loc_pass=FALSE
      }
      
      return(c(CpG=cpg, nmod=nrow(modes), loc_pass=loc_pass, 
               loc0=loc012[1], loc1=loc012[2], loc2=loc012[3], 
               dens0=dens012[1], dens1=dens012[2], dens2=dens012[3]))
    }, error = function(e) return(paste0("Escape ", cpg, " with error: ", e)))
  }
  stopImplicitCluster()
  stopCluster(cl)
  modRes <- as.data.frame(do.call(rbind, modRes))
  modRes$nmod <- as.numeric(as.character(modRes$nmod))
  modRes$loc0 <- as.numeric(as.character(modRes$loc0))
  modRes$loc1 <- as.numeric(as.character(modRes$loc1))
  modRes$loc2 <- as.numeric(as.character(modRes$loc2))
  modRes$dens0 <- as.numeric(as.character(modRes$dens0))
  modRes$dens1 <- as.numeric(as.character(modRes$dens1))
  modRes$dens2 <- as.numeric(as.character(modRes$dens2))
  rownames(modRes) <- modRes$CpG
  modRes
}
