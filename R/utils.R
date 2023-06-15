.onAttach <- function(libname, pkgname) {
  message("Loading package: ", pkgname)
  message("Located in library: ", libname)
  # packageStartupMessage("Welcome")
  
  suppressMessages({
    library(MethylGenotyper)
    library(minfi)
    library(tidyverse)
    library(foreach)
    library(doParallel)
    library(HardyWeinberg)
    library(multimode)
    library(rlist)
    library(stats4)
    library(ggplot2)
    library(ggpubr)
    library(ggsci)
  })
}
