
#' Call genotypes based on EM algorithm and Bayesian approach
#' 
#' The Expectation–maximization (EM) algorithm is used to fit a mixture of three beta distributions representing the three genotypes (AA, AB, and BB) and one uniform distribution representing the outliers (adapted from ewastools). Then, the Bayesian approach is used to get the genotype probabilities, with AFs of the matched population (the 1000 Genomes Project, 1KGP) being used to infer priors.
#'
#' @param RAI A matrix of RAI (Ratio of Alternative allele Intensity) for probes. Provide probes as rows and samples as columns.
#' @param pop Population to be used to extract AFs. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @param maxiter Maximal number of iterations for the EM algorithm.
#' @return  A list containing
#' \item{RAI}{Ratio of Alternative allele Intensity}
#' \item{shapes}{Shapes of the mixed beta distributions}
#' \item{iterations}{Iterations}
#' \item{GP}{Posterior probabilities for the three genotypes}
#' @export
call_genotypes_bayesian <- function(RAI, pop, type, maxiter=50){
  print(paste(Sys.time(), "Running EM to fit beta distributions for RAI values."))
  # Initialize
  assignments <- RAI %>% as.data.frame %>%
    mutate(Probe=rownames(.)) %>%
    tidyr::gather(key="Sample", value="RAI", -Probe) %>% tibble
  N <- nrow(assignments)
  
  # EM
  priors <- rep(1/3, 3)
  shape1_1 <- 5; shape2_1 <- 60
  shape1_2 <- 30; shape2_2 <- 30
  shape1_3 <- 60; shape2_3 <- 5
  U <- 0.01
  outliers <- rep(U, N)
  GP <- NA

  e_step <- function(){
    GP <- (1 - U) * cbind(
      priors[1] * dbeta(assignments$RAI, shape1_1, shape2_1, log=F),
      priors[2] * dbeta(assignments$RAI, shape1_2, shape2_2, log=F),
      priors[3] * dbeta(assignments$RAI, shape1_3, shape2_3, log=F)
    )
    tmp <- rowSums(GP)
    GP <<- GP / tmp
    outliers <<- U / (U + tmp)
    logLik <- sum(log(U + tmp))
    return(logLik)
  }
  
  m_step <- function(){
    GP <- GP * (1 - outliers)
    priors <- colSums(GP)
    priors <<- priors / sum(priors)
    U <<- sum(outliers) / N
    # Moments estimator
    s1 <- eBeta(assignments$RAI, GP[,1])
    s2 <- eBeta(assignments$RAI, GP[,2])
    s3 <- eBeta(assignments$RAI, GP[,3])
    shape1_1 <<- s1$shape1; shape1_2 <<- s2$shape1; shape1_3 <<- s3$shape1
    shape2_1 <<- s1$shape2; shape2_2 <<- s2$shape2; shape2_3 <<- s3$shape2
    return(NA)
  }
  
  iterations <- list()
  gain <- Inf
  i <- 1
  while(i < maxiter & gain > 1e-4){
    logLik <- e_step()
    m_step()
    shapes <- as.data.frame(matrix(
      c(shape1_1, shape1_2, shape1_3, shape2_1, shape2_2, shape2_3),
      nrow=3,
      dimnames=list(paste0("Cluster", 1:3), c("shape1", "shape2"))
    ))
    iterations[[i]] <- list(
      priors = priors,
      shapes = shapes,
      U = U,
      logLik = logLik
    )
    print(paste0("EM: Iteration ", i, ", Minus log-likelihood: ", round(logLik, 6)))
    # print("Prior probabilities of the three genotypes: ")
    # print(priors)
    # print("Shapes:")
    # print(shapes)
    # print("Outlier:")
    # print(U)
    if(i>1){
      gain = logLik - iterations[[i-1]]$logLik
    }
    i=i+1
  }
  finalClusters <- iterations[[i-1]]
  
  # Get posterior genotype probabilities
  print(paste(Sys.time(), "Running the Bayesian approach to get posterior genotype probabilities."))
  probe2af <- get_AF(pop, type)
  AF <- matrix(rep(probe2af[rownames(RAI)], ncol(RAI)), ncol=ncol(RAI), dimnames=list(rownames(RAI), colnames(RAI)))
  GP <- get_GP(RAI, finalClusters$shapes[, c("shape1", "shape2")], bayesian=TRUE, AF)
  
  # return
  list(RAI=RAI, 
       shapes=finalClusters$shapes, 
       priors=finalClusters$priors, 
       outlier=finalClusters$U, 
       iterations=iterations, 
       GP=GP)
}

#' Moments estimator for beta distribution (adapted from ewastools)
#' @export
eBeta = function(x,w){
  n = length(w)
  w = n*w/sum(w)
  sample.mean =  mean(w*x)
  sample.var  = (mean(w*x^2)-sample.mean^2) * n/(n-1)
  v = sample.mean * (1-sample.mean)
  if (sample.var < v){
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
#' @param pop Population to be used to extract AFs. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @return A vector of AFs
#' @export
get_AF <- function(pop="EAS", type){
  if(!(pop %in% c("EAS", "AMR", "AFR", "EUR", "SAS", "ALL"))){
    print("ERROR: get_AF: Wrong pop value supplied!"); return(NA)
  }
  if(type=="snp_probe"){
    data(probeInfo_snp)
    probe2af <- probeInfo_snp[, paste0(pop, "_AF")]
    names(probe2af) <- probeInfo_snp$CpG
  }else if(type=="typeI_ccs_probe"){
    data(probeInfo_typeI)
    probe2af <- probeInfo_typeI[, paste0(pop, "_AF")]
    names(probe2af) <- probeInfo_typeI$CpG
  }else if(type=="typeII_ccs_probe"){
    data(probeInfo_typeII)
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
#' @param RAI A MxN matrix of RAI (Ratio of Alternative allele Intensity). Provide probes as rows and samples as columns.
#' @param shapes A data frame (3x2) containing the two shapes for beta distributions of the three clusters. 
#' @param bayesian Use the Bayesian approach or not.
#' @param AF A MxN matrix of AFs. Provide SNPs as rows and samples as columns.
#' @return  A list containing
#' \item{pAA}{Posterior genotype probability of AA}
#' \item{pAB}{Posterior genotype probability of AB}
#' \item{pBB}{Posterior genotype probability of BB}
#' @export
get_GP <- function(RAI, shapes, bayesian=TRUE, AF){
  probes <- intersect(rownames(RAI), rownames(AF))
  samples <- intersect(colnames(RAI), colnames(AF))
  RAI <- RAI[probes, samples]
  shapes$mean <- shapes$shape1 / (shapes$shape1 + shapes$shape2)
  shapes <- as.matrix(shapes[order(shapes$mean),]) # row1 to row3: AA, AB, BB
  pD_AA <- apply(RAI, 1:2, function(x) dbeta(x, shapes[1, 1], shapes[1, 2])) # probability of data given genotype AA
  pD_AB <- apply(RAI, 1:2, function(x) dbeta(x, shapes[2, 1], shapes[2, 2]))
  pD_BB <- apply(RAI, 1:2, function(x) dbeta(x, shapes[3, 1], shapes[3, 2]))
  if(bayesian==TRUE){
    AF <- AF[probes, samples]
    pAA_prior <- (1 - AF) ^2 # Prior genotype probability of AA
    pAB_prior <- 2 * AF * (1 - AF)
    pBB_prior <- AF ^2
    pD <- pD_AA * pAA_prior + pD_AB * pAB_prior + pD_BB * pBB_prior
    pAA <- pD_AA * pAA_prior / pD
    pAB <- pD_AB * pAB_prior / pD
    pBB <- pD_BB * pBB_prior / pD
  }else{
    pD <- pD_AA + pD_AB + pD_BB
    pAA <- pD_AA / pD
    pAB <- pD_AB / pD
    pBB <- pD_BB / pD
  }
  return(list(pAA=pAA, pAB=pAB, pBB=pBB))
}


