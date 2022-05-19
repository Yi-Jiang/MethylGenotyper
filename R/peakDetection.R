
#' Filter probes according to the number and position of peaks
#' 
#' @param AB_geno Matrix of allelic balances.
#' @param cpu Number of CPU.
#' @return Probe filtering results
#' @export
getMod <- function(AB_geno, cpu=1){
  # Estimate number of modes for each probe
  n <- apply(AB_geno, 1, function(x) (nmodes(x, 0.05, lowsup=0, uppsup=1)))
  mod <- tibble(Name=names(n), nmod=n)
  mod$nmod2 <- mod$nmod
  mod[mod$nmod2>3, "nmod2"] <- 3
  
  # Estimate mode location for each probe
  cl<- makeCluster(cpu)
  registerDoParallel(cl) 
  out <- foreach(i=mod$Name, .packages=c("tidyverse","multimode")) %dopar% {
    tryCatch({
      a <- locmodes(AB_geno[i, ], mod0=filter(mod, Name==i)$nmod2)
      a
    }, error = function(e) return(paste0("The variable '", i, "'", " caused the error: '", e, "'")))
  }
  stopImplicitCluster()
  stopCluster(cl)
  names(out) <- mod$Name
  out2 <- out[filter(mod, nmod2==2)$Name]
  out3 <- out[filter(mod, nmod2==3)$Name]
  
  # Filter mode density height (>0.1)----
  height2 <- as.data.frame(do.call(rbind, list.map(out2, fvalue))) %>% select(V1, V3)
  height3 <- as.data.frame(do.call(rbind, list.map(out3, fvalue))) %>% select(V1, V3, V5)
  height <- bind_rows(height2, height3)
  dat <- gather(height, `V1`,`V3`,`V5`, key="mode", value="h")
  hfilter <- rowSums(height>0.1, na.rm=T)>1
  nhfilter <- rowSums(height>0.1, na.rm=T)
  mod <- left_join(mod, tibble(Name=names(hfilter), h_0.1=hfilter, nh_0.1=nhfilter))
  
  # Filter mode location ----
  lo2 <- as.data.frame(do.call(rbind, list.map(out2, locations))) %>% select(V1, V3)
  lo3 <- as.data.frame(do.call(rbind, list.map(out3, locations))) %>% select(V1, V3, V5)
  k.lo2 <- filter(lo2, (V1<0.3&V3>0.3&V3<0.7)|(V1>0.3&V1<0.7&V3>0.7)) 
  k.lo3 <- filter(lo3, V1<0.3, V3>0.3, V3<0.7, V5>0.7) 
  k.lo <- bind_rows(k.lo2, k.lo3)
  mod <- mutate(mod, loc_pass=(Name %in% rownames(k.lo)))
  
  mod
}



