
#' Get example IDAT file list
#' 
#' @param type. One of "EPIC" and "450K"
#' @return A data frame of the IDAT file list
get_target <- function(type="EPIC"){
  if(type=="EPIC"){
    library(minfiDataEPIC)
    baseDir <- system.file("extdata", package = "minfiDataEPIC")
  }else if(type=="450K" | type=="450k"){
    library(minfiData)
    baseDir <- system.file("extdata", package = "minfiData")
  }else{
    print("Please specify one of EPIC and 450K.")
    return
  }
  target <- list()
  filelist <- list.files(baseDir)
  filelist <- filelist[-length(filelist)]
  for(i in filelist){
    path <- file.path(baseDir, i)
    prefix <- sub("_Grn.idat", "", list.files(path)[1])
    target[[i]] <- c(Sample_Name=i, Basename=paste(path, prefix, sep="/"))
  }
  target <- as.data.frame(do.call(rbind, target))
  rownames(target) <- 1:nrow(target)
  target
}
