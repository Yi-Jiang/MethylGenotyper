
#' Get example IDAT file list
#' 
#' @param platform One of "EPIC" and "450K"
#' @return A data frame of the IDAT file list
#' @export
get_target <- function(platform="EPIC"){
  if(grepl("EPIC", platform, ignore.case=TRUE)){
    suppressMessages(library(minfiDataEPIC))
    baseDir <- system.file("extdata", package = "minfiDataEPIC")
  }else if(grepl("450", platform, ignore.case=TRUE)){
    suppressMessages(library(minfiData))
    baseDir <- system.file("extdata", package = "minfiData")
  }else{
    print("Please specify the platform type: EPIC or 450K.")
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
