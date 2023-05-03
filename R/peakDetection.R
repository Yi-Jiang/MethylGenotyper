
#' Filter probes according to the number and position of peaks
#' 
#' @param RAI Matrix of RAI (Ratio of Alternative allele Intensity).
#' @param cpu Number of CPU.
#' @return Probe filtering results
#' @export
getMod <- function(RAI, cpu=1){
  print(paste(Sys.time(), "Running mod test."))
  # Estimate number of modes for each probe
  n <- apply(RAI, 1, function(x) (nmodes(x, 0.05, lowsup=0, uppsup=1)))
  mod <- tibble(Name=names(n), nmod=n)
  mod$nmod2 <- mod$nmod
  mod[mod$nmod2>3, "nmod2"] <- 3
  
  # Estimate mode location for each probe
  cl<- makeCluster(cpu)
  registerDoParallel(cl) 
  out <- foreach(i=mod$Name, .packages=c("tidyverse","multimode")) %dopar% {
    tryCatch({
      a <- locmodes(RAI[i, ], mod0=dplyr::filter(mod, Name==i)$nmod2)
      a
    }, error = function(e) return(paste0("The variable '", i, "'", " caused the error: '", e, "'")))
  }
  stopImplicitCluster()
  stopCluster(cl)
  names(out) <- mod$Name
  out2 <- out[dplyr::filter(mod, nmod2==2)$Name]
  out3 <- out[dplyr::filter(mod, nmod2==3)$Name]
  
  # Filter mode density height (>0.1)
  if(length(out2)>0){
    height2 <- as.data.frame(do.call(rbind, list.map(out2, fvalue))) %>% dplyr::select(V1, V3)
  }else{
    height2 <- c()
  }
  if(length(out3)>0){
    height3 <- as.data.frame(do.call(rbind, list.map(out3, fvalue))) %>% dplyr::select(V1, V3, V5)
  }else{
    height3 <- c()
  }
  height <- bind_rows(height2, height3)
  hfilter <- rowSums(height>0.1, na.rm=T)>1
  nhfilter <- rowSums(height>0.1, na.rm=T)
  mod <- left_join(mod, tibble(Name=names(hfilter), h_0.1=hfilter, nh_0.1=nhfilter))
  
  # Filter mode location
  if(length(out2)>0){
    lo2 <- as.data.frame(do.call(rbind, list.map(out2, locations))) %>% dplyr::select(V1, V3)
    k.lo2.01 <- data.frame(dplyr::filter(lo2, V1<0.3&V3>0.3&V3<0.7), V5=NA_real_)
    k.lo2.12 <- data.frame(V0=NA_real_, dplyr::filter(lo2, V1>0.3&V1<0.7&V3>0.7)); colnames(k.lo2.12) <- c("V1", "V3", "V5")
  }else{
    k.lo2.01 <- c()
    k.lo2.12 <- c()
  }
  if(length(out3)>0){
    lo3 <- as.data.frame(do.call(rbind, list.map(out3, locations))) %>% dplyr::select(V1, V3, V5)
    k.lo3 <- dplyr::filter(lo3, V1<0.3, V3>0.3, V3<0.7, V5>0.7)
  }else{
    k.lo3 <- c()
  }
  k.lo <- bind_rows(k.lo2.01, k.lo2.12, k.lo3)
  colnames(k.lo) <- c("loc0", "loc1", "loc2")
  k.lo <- data.frame(Name=rownames(k.lo), k.lo)
  mod <- mutate(mod, loc_pass=(Name %in% rownames(k.lo))) %>% left_join(k.lo, by="Name")
  mod <- as.data.frame(mod)
  rownames(mod) <- mod$Name

  mod
}



