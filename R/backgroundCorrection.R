
#' Noob and dye-bias correction
#' 
#' @param target A data frame of two columns: Sample_Name, Basename, where Basename tells the location of IDAT files.
#' @param mnfst Manifest file. Required columns: Name, AddressA_ID (numeric), AddressB_ID (numeric), Infinium_Design_Type, and Color_Channel
#' @param cpu Number of CPU.
#' @return A list of noob and dye-bias corrected signals containing:
#' \item{AR}{ - A matrix of probeA signals in Red channel}
#' \item{AG}{ - A matrix of probeA signals in Green channel}
#' \item{BR}{ - A matrix of probeB signals in Red channel}
#' \item{BG}{ - A matrix of probeB signals in Green channel}
#' @export
correct_noob_dye <- function(target, mnfst, cpu=1){
    data(probelist)
    cl <- makeCluster(cpu)
    registerDoParallel(cl)
    rgData_list <- foreach (sp=c(1:nrow(target)), .packages="tidyverse", .export=c("probelist", "backgroundCorrectionNoobFit", "normExpSignal")) %dopar% {
        rgSet = minfi::read.metharray.exp(targets=target[sp,])
        green <- minfi::getGreen(rgSet)[,1] # vector: IBG IAG IIG
        red <- minfi::getRed(rgSet)[,1]   # vector: IBR IAR IIR
        green[green == 0] <- 1
        red[red == 0] <- 1
        
        ## format in sdg
        rgData0 <- left_join(mnfst, tibble(AG=green, AddressA_ID=as.integer(names(green))), by=c("AddressA_ID")) %>% 
            left_join(tibble(AR=red, AddressA_ID=as.integer(names(red))), by=c("AddressA_ID")) %>% 
            left_join(tibble(BG=green, AddressB_ID=as.integer(names(green))), by=c("AddressB_ID")) %>% 
            left_join(tibble(BR=red, AddressB_ID=as.integer(names(red))), by=c("AddressB_ID"))
        dG <- dplyr::filter(rgData0, Color_Channel=="Grn")
        dR <- dplyr::filter(rgData0, Color_Channel=="Red")
        InfII <- dplyr::filter(rgData0, Infinium_Design_Type=="II")
        
        ## foreground (in band, ib)
        ibG <- c(dG$BG, dG$AG, InfII$AG)
        ibR <- c(dR$BR, dR$AR, InfII$AR)
        ibG[ibG == 0] <- 1
        ibR[ibR == 0] <- 1
        
        ## background (out of band, oob)
        oobG <- c(dR$BG, dR$AG)
        oobR <- c(dG$BR, dG$AR)
        oobR[oobR == 0] <- 1
        oobG[oobG == 0] <- 1
        
        ## Noob for Green channel
        fitG <- backgroundCorrectionNoobFit(ibG, oobG)
        rgData0$BG_noob <- normExpSignal(fitG$mu, fitG$sigma, fitG$alpha, rgData0$BG) + 15
        rgData0$AG_noob <- normExpSignal(fitG$mu, fitG$sigma, fitG$alpha, rgData0$AG) + 15
        
        ## Noob for red channel
        fitR <- backgroundCorrectionNoobFit(ibR, oobR)
        rgData0$BR_noob <- normExpSignal(fitR$mu, fitR$sigma, fitR$alpha, rgData0$BR) + 15
        rgData0$AR_noob <- normExpSignal(fitR$mu, fitR$sigma, fitR$alpha, rgData0$AR) + 15
        
        ## Only keep candidate SNP probes, Type I CCS probes, and Type II CCS probes
        rgData0 <- rgData0[rgData0$Name %in% probelist$CpG,]
        
        ## Fit linear regression of red and green channel
        control_probes <- minfi::getProbeInfo(rgSet, type = "Control")
        controls <- as.data.frame(control_probes) %>% 
            left_join(tibble(G=green, Address=names(green))) %>% 
            left_join(tibble(R=red, Address=names(red)))
        Ai = dplyr::filter(controls, Type=="NORM_A") %>% arrange(ExtendedType) # 27
        Ti = dplyr::filter(controls, Type=="NORM_T") %>% arrange(ExtendedType) # 58
        Gi = dplyr::filter(controls, Type=="NORM_G") %>% arrange(ExtendedType) # 27
        Ci = dplyr::filter(controls, Type=="NORM_C") %>% arrange(ExtendedType) # 58
        x <- log(c(Gi$G, Ci$G)) # Green
        y <- log(c(Ai$R, Ti$R)) # Red
        keep = !is.na(y) & !is.na(x) & is.finite(x) & is.finite(y)
        x = x[keep]
        y = y[keep]
        m = mblm::mblm(y~x, repeated=FALSE)
        
        ## dye bias correction for green channel
        dG.noob <- dplyr::filter(rgData0, Color_Channel=="Grn")
        dR.noob <- dplyr::filter(rgData0, Color_Channel=="Red")
        InfII.noob <- dplyr::filter(rgData0, Infinium_Design_Type=="II")
        dG.noob$BG_noob_dye = exp(coef(m)[1] + log(dG.noob$BG_noob) * coef(m)[2])
        dG.noob$AG_noob_dye = exp(coef(m)[1] + log(dG.noob$AG_noob) * coef(m)[2])
        dR.noob$BG_noob_dye = exp(coef(m)[1] + log(dR.noob$BG_noob) * coef(m)[2])  # oob
        dR.noob$AG_noob_dye = exp(coef(m)[1] + log(dR.noob$AG_noob) * coef(m)[2])  # oob
        InfII.noob$BG_noob_dye = NA
        InfII.noob$AG_noob_dye = exp(coef(m)[1] + log(InfII.noob$AG_noob) * coef(m)[2])
        rgData0 <- rbind(dG.noob, dR.noob, InfII.noob)
        rgData0
    }
    stopImplicitCluster()
    stopCluster(cl)
    
    ## format results
    rgData <- list("AG"=list(), "AR"=list(), "BG"=list(), "BR"=list())
    sampleNames <- target$Sample_Name
    for(i in 1:length(rgData_list)){
        sam <- sampleNames[i]
        rgData[["AG"]][[sam]] <- rgData_list[[i]]$AG_noob_dye
        rgData[["AR"]][[sam]] <- rgData_list[[i]]$AR_noob
        rgData[["BG"]][[sam]] <- rgData_list[[i]]$BG_noob_dye
        rgData[["BR"]][[sam]] <- rgData_list[[i]]$BR_noob
    }
    rgData[["AG"]] <- do.call(cbind, rgData[["AG"]]); rownames(rgData[["AG"]]) <- rgData_list[[1]]$Name
    rgData[["AR"]] <- do.call(cbind, rgData[["AR"]]); rownames(rgData[["AR"]]) <- rgData_list[[1]]$Name
    rgData[["BG"]] <- do.call(cbind, rgData[["BG"]]); rownames(rgData[["BG"]]) <- rgData_list[[1]]$Name
    rgData[["BR"]] <- do.call(cbind, rgData[["BR"]]); rownames(rgData[["BR"]]) <- rgData_list[[1]]$Name
    rgData
}

#' Fit Normal and exponential distributions (adapted from `SeSAMe`)
#'
#' @param ib Foreground signals.
#' @param bg Background signals.
#' @return mu and sigma by fitting background signals with normal distribution.
#' @return alpha by fitting foreground signals with exponential distribution.
#' @export
backgroundCorrectionNoobFit <- function(ib, bg){
    e <- MASS::huber(bg)
    mu <- e$mu
    sigma <- e$s
    alpha <- pmax(MASS::huber(ib)$mu-mu, 10)
    list(mu = mu, sigma = sigma, alpha = alpha)
}

#' Normal-exponential deconvolution (adapted from `SeSAMe`)
#' 
#' @param mu,sigma Background signal parameters returned by `backgroundCorrectionNoobFit`.
#' @param alpha Foreground signal parameters returned by `backgroundCorrectionNoobFit`.
#' @param x Foreground signals to be corrected.
#' @return The conditional expectation of the signal given the observed foreground and background.
#' @export
normExpSignal <- function(mu, sigma, alpha, x){
    sigma2 <- sigma * sigma
    if (alpha <= 0)
        stop("alpha must be positive")
    if (sigma <= 0)
        stop("sigma must be positive")
    mu.sf <- x - mu - sigma2/alpha
    signal <- mu.sf + sigma2 * exp(
        dnorm(0, mean = mu.sf, sd = sigma, log = TRUE) -
            pnorm(
                0, mean = mu.sf, sd = sigma,
                lower.tail = FALSE, log.p = TRUE))
    o <- !is.na(signal)
    if (any(signal[o] < 0)) {
        warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
        signal[o] <- pmax(signal[o], 1e-06)
    }
    signal
}

