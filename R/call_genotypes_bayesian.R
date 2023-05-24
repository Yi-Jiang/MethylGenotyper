
#' Call genotypes based on EM algorithm and Bayesian approach
#' 
#' The Expectation–maximization (EM) algorithm is used to fit a mixture of three beta distributions representing the three genotypes (AA, AB, and BB) and one uniform distribution representing the outliers (adapted from ewastools). Then, the Bayesian approach is used to get the genotype probabilities, with AFs of the matched population (the 1000 Genomes Project, 1KGP) being used to infer priors.
#'
#' @param RAI A matrix of RAI (Ratio of Alternative allele Intensity) for probes. Provide probes as rows and samples as columns.
#' @param pop Population to be used to extract AFs. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @param maxiter Maximal number of iterations for the EM algorithm.
#' @param bayesian Use the Bayesian approach to calculate posterior genotype probabilities.
#' @param platform EPIC or 450K.
#' @param verbose Verbose mode: 0/1/2.
#' @return  A list containing
#' \item{RAI}{Ratio of Alternative allele Intensity}
#' \item{shapes}{Shapes of the mixed beta distributions}
#' \item{weights}{Prior probabilities that the RAI values belong to one of the three genotypes}
#' \item{U}{Overall probability of RAI values being outlier}
#' \item{outliers}{Probability of each RAI value being outlier}
#' \item{logLik}{Log-likelihood}
#' \item{GP}{Genotype probabilities of the three genotypes}
#' @export
call_genotypes_bayesian <- function(RAI, pop, type, maxiter=50, bayesian=TRUE, platform="EPIC", verbose=1){
  # Fit mixed beta distribution based on EM
  print(paste(Sys.time(), "Running EM to fit beta distributions for RAI values."))
  finalClusters <- fit_beta_em(RAI, maxiter=maxiter, verbose=verbose)
  GP <- finalClusters$GP

  # Get posterior genotype probabilities using a Bayesian approach
  if(bayesian){
    print(paste(Sys.time(), "Calculating genotype probabilities using a Bayesian approach."))
    probe2af <- get_AF(pop=pop, type=type, platform=platform)
    AF <- matrix(rep(probe2af[rownames(RAI)], ncol(RAI)), ncol=ncol(RAI), dimnames=list(rownames(RAI), colnames(RAI)))
    GP <- get_GP_bayesian(GP$pAA, GP$pAB, GP$pBB, AF)
  }
  
  # return
  list(RAI=RAI, 
       shapes = finalClusters$shapes, 
       weights = finalClusters$weights, 
       U = finalClusters$U, 
       outliers = finalClusters$outliers,
       logLiks = finalClusters$logLiks,
       GP = GP)
}

#' Estimate mixed beta distribution parameters based on EM algorithm
#' 
#' The Expectation–maximization (EM) algorithm is used to fit a mixture of three beta distributions representing the three genotypes (AA, AB, and BB) and one uniform distribution representing the outliers (adapted from ewastools).
#'
#' @param RAI A matrix of RAI (Ratio of Alternative allele Intensity) for probes. Provide probes as rows and samples as columns.
#' @param maxiter Maximal number of iterations for the EM algorithm.
#' @param verbose Verbose mode: 0/1/2.
#' @return  A list containing
#' \item{shapes}{Shapes of the mixed beta distributions}
#' \item{weights}{Prior probabilities that the RAI values belong to one of the three genotypes}
#' \item{U}{Overall probability of RAI values being outlier}
#' \item{outliers}{Probability of each RAI value being outlier}
#' \item{logLik}{Log-likelihood}
#' \item{GP}{Genotype probabilities of the three genotypes}
#' @export
fit_beta_em <- function(RAI, maxiter=50, verbose=1){
  # Initialize
  RAI_flat <- as.vector(RAI) # 1-dimension
  m <- nrow(RAI) # number of probes
  n <- ncol(RAI) # number of samples
  N <- length(RAI_flat)
  
  # EM
  weights <- rep(1/3, 3)
  shapes <- as.data.frame(matrix(c(5, 30, 60, 60, 30, 5), nrow=3, dimnames=list(paste0("Cluster", 0:2), c("shape1", "shape2"))))
  U <- 0.01
  outliers <- rep(U, N)

  e_step <- function(){
    GP <- (1 - U) * cbind(
      weights[1] * dbeta(RAI_flat, shapes[1,1], shapes[1,2], log=F),
      weights[2] * dbeta(RAI_flat, shapes[2,1], shapes[2,2], log=F),
      weights[3] * dbeta(RAI_flat, shapes[3,1], shapes[3,2], log=F)
    )
    tmp <- rowSums(GP, na.rm=T)
    GP <<- GP / tmp # If RAI=0，the "dbeta" function above will produce three zeros, and leads to NA here.
    outliers <<- U / (U + tmp)
    logLik <- sum(log(U + tmp))
    return(logLik)
  }
  
  m_step <- function(){
    GP <- GP * (1 - outliers)
    weights <- colSums(GP)
    weights <<- weights / sum(weights, na.rm=T)
    U <<- sum(outliers) / N
    # Moments estimator
    s1 <- eBeta(RAI_flat, GP[,1])
    s2 <- eBeta(RAI_flat, GP[,2])
    s3 <- eBeta(RAI_flat, GP[,3])
    shapes[1,1] <<- s1$shape1; shapes[2,1] <<- s2$shape1; shapes[3,1] <<- s3$shape1
    shapes[1,2] <<- s1$shape2; shapes[2,2] <<- s2$shape2; shapes[3,2] <<- s3$shape2
    return(NA)
  }
  
  logLiks <- c()
  gain <- Inf
  i <- 1
  while(i <= maxiter & gain > 1e-4){
    logLik <- e_step()
    m_step()
    logLiks <- c(logLiks, logLik)
    if(verbose>=1){
      print(paste0("EM: Iteration ", i, ", log-likelihood: ", round(logLik, 6)))
    }
    if(i>1){
      gain = logLik - logLiks[i-1]
    }
    i=i+1
  }
  names(logLiks) <- paste0("Iter", 1:length(logLiks))
  outliers <- matrix(outliers, nrow=m, ncol=n, dimnames=list(rownames(RAI), colnames(RAI)))
  GP_AA <- matrix(GP[,1], nrow=m, ncol=n, dimnames=list(rownames(RAI), colnames(RAI)))
  GP_AB <- matrix(GP[,2], nrow=m, ncol=n, dimnames=list(rownames(RAI), colnames(RAI)))
  GP_BB <- matrix(GP[,3], nrow=m, ncol=n, dimnames=list(rownames(RAI), colnames(RAI)))
  iteration <- list(
    weights = weights,
    shapes = shapes,
    U = U,
    outliers = outliers,
    logLiks = logLiks,
    GP = list(pAA = GP_AA, pAB = GP_AB, pBB = GP_BB)
  )
  iteration
}

#' Moments estimator for beta distribution (adapted from ewastools)
#' @export
eBeta = function(x,w){
  n = length(w)
  w = n*w/sum(w)
  sample.mean =  mean(w*x)
  sample.var  = (mean(w*x^2)-sample.mean^2) * n/(n-1)
  v = sample.mean * (1-sample.mean)
  if (is.na(sample.var) | is.na(v)) {
    shape1 = Inf
    shape2 = Inf
  } else if (sample.var < v){
    shape1 = sample.mean * (v/sample.var - 1)
    shape2 = (1 - sample.mean) * (v/sample.var - 1)
  } else {
    shape2 = sample.mean * (v/sample.var - 1)
    shape1 = (1 - sample.mean) * (v/sample.var - 1)
  }
  list(shape1 = shape1, shape2 = shape2)
}

#' Extract AFs from matching population in the 1000 Genomes Project (1KGP)
#'
#' @param platform EPIC or 450K.
#' @param pop Population to be used to extract AFs. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @param platform EPIC or 450K.
#' @return A vector of AFs
#' @export
get_AF <- function(pop="EAS", type, platform="EPIC"){
  if(!(pop %in% c("EAS", "AMR", "AFR", "EUR", "SAS", "ALL"))){
    print("ERROR: get_AF: Wrong pop value supplied!"); return(NA)
  }
  if(type=="snp_probe"){
    if(platform=="EPIC"){
      data(probeInfo_snp)
    }else{
      data(probeInfo_snp_450K); probeInfo_snp <- probeInfo_snp_450K
    }
    probe2af <- probeInfo_snp[, paste0(pop, "_AF")]
    names(probe2af) <- probeInfo_snp$CpG
  }else if(type=="typeI_ccs_probe"){
    if(platform=="EPIC"){
      data(probeInfo_typeI)
    }else{
      data(probeInfo_typeI_450K); probeInfo_typeI <- probeInfo_typeI_450K
    }
    probe2af <- probeInfo_typeI[, paste0(pop, "_AF")]
    names(probe2af) <- probeInfo_typeI$CpG
  }else if(type=="typeII_ccs_probe"){
    if(platform=="EPIC"){
      data(probeInfo_typeII)
    }else{
      data(probeInfo_typeII_450K); probeInfo_typeII <- probeInfo_typeII_450K
    }
    probe2af <- probeInfo_typeII[, paste0(pop, "_AF")]
    names(probe2af) <- probeInfo_typeII$CpG
  }else{
    print("ERROR: get_AF: Wrong type value supplied!"); return(NA)
  }
  probe2af
}

#' Infer posterior genotype probabilities based on the Bayesian approach
#'
#' Prior genotype probabilities were inferred from AFs. The AFs can be in population level or individual-specific level. For population level AFs, they can be extracted from the matched population in the 1000 Genomes Project (1KGP). For individual-specific AFs, they can be calculated according to the top four PCs.
#'
#' @param pD_AA A MxN matrix of AA genotype probabilities. Provide SNPs as rows and samples as columns.
#' @param pD_AB A MxN matrix of AB genotype probabilities. Provide SNPs as rows and samples as columns.
#' @param pD_BB A MxN matrix of BB genotype probabilities. Provide SNPs as rows and samples as columns.
#' @param AF A MxN matrix of AFs. Provide SNPs as rows and samples as columns.
#' @return  A list containing
#' \item{pAA}{Posterior genotype probability of AA}
#' \item{pAB}{Posterior genotype probability of AB}
#' \item{pBB}{Posterior genotype probability of BB}
#' @export
get_GP_bayesian <- function(pD_AA, pD_AB, pD_BB, AF){
  probes <- intersect(rownames(pD_AA), rownames(AF))
  samples <- intersect(colnames(pD_AA), colnames(AF))
  pD_AA <- pD_AA[probes, samples]
  pD_AB <- pD_AB[probes, samples]
  pD_BB <- pD_BB[probes, samples]
  AF <- AF[probes, samples]
  pAA_prior <- (1 - AF) ^2 # Prior genotype probability of AA
  pAB_prior <- 2 * AF * (1 - AF)
  pBB_prior <- AF ^2
  pD <- pD_AA * pAA_prior + pD_AB * pAB_prior + pD_BB * pBB_prior
  pAA <- pD_AA * pAA_prior / pD
  pAB <- pD_AB * pAB_prior / pD
  pBB <- pD_BB * pBB_prior / pD
  return(list(pAA=pAA, pAB=pAB, pBB=pBB))
}

