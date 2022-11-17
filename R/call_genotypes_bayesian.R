
#' Fit beta distribution based on Maximum likelihood estimation
#' 
#' @param x A vector of RAI values.
#' @return A data frame with one row and four columns:
#' \item{shape1}{The first parameter for beta distribution}
#' \item{shape2}{The second parameter for beta distribution}
#' \item{mLL}{Minus log-likelihood scaled by number of data points}
#' \item{number}{The length of x}
#' @export
fit_beta_mle <- function(x){
  minuslogL <- function(shape1, shape2){-sum(dbeta(x, shape1, shape2, log=T))}
  m <- stats4::mle(minuslogL, start=list(shape1=3, shape2=3), method="L-BFGS-B", lower=c(0.001, 0.001))
  shape <- stats4::coef(m)
  params <- tibble(shape1=shape[1], shape2=shape[2], mLL=m@min/length(x), number=length(x))
  params
}

#' Get clusters based on the Expectation–maximization (EM) algorithm
#' 
#' @param assignments A data frame containing:
#' \item{Probe}{Probe ID}
#' \item{Sample}{Sample ID}
#' \item{RAI}{Ratio of Alternative allele Intensity}
#' \item{Cluster}{Randomly assigned three clusters indicating the three genotypes}
#' @return A list of two elements:
#' \item{assignments}{A data frame with re-assigned clusters}
#' \item{fits}{A data frame containing the two shape parameters for the three re-assigned clusters}
#' @export
iterate_em <- function(assignments){
  fits <- assignments %>%
    group_by(Cluster) %>%
    do(mutate(fit_beta_mle(.$RAI))) %>%
    ungroup() %>%
    mutate(prior = number / sum(number))
  
  assignments <- assignments %>%
    dplyr::select(Probe:RAI) %>%
    crossing(fits) %>%
    mutate(likelihood = prior * dbeta(RAI, shape1, shape2)) %>%
    group_by(Probe, Sample) %>%
    top_n(1, likelihood) %>%
    ungroup()
  
  list(assignments = assignments, fits = fits)
}

#' Extract AFs from matching population in the 1000 Genomes Project (1KGP)
#'
#' @param RAI A matrix of RAI (Ratio of Alternative allele Intensity) for probes. Provide probes as rows and samples as columns.
#' @param pop Population to be used to extract AFs. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @return A vector of AFs
#' @export
get_AF <- function(RAI, pop="EAS", type){
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
    print("ERROR: get_GP: Wrong pop value supplied!"); return(NA)
  }
  probe2af <- probe2af[rownames(RAI)]
  probe2af
}

#' Infer posterior genotype probabilities based on the Bayesian approach
#'
#' Prior genotype probabilities were inferred from AFs of matching population in the 1000 Genomes Project (1KGP).
#'
#' @param RAI A matrix of RAI (Ratio of Alternative allele Intensity) for probes. Provide probes as rows and samples as columns.
#' @param shapes A data frame (3x2) containing the two shapes for beta distributions of the three clusters. 
#' @param probe2af A vector of AFs extracted from matching population in 1KGP.
#' @return  A list containing
#' \item{pAA}{Posterior genotype probability of AA}
#' \item{pAB}{Posterior genotype probability of AB}
#' \item{pBB}{Posterior genotype probability of BB}
#' @export
get_GP <- function(RAI, shapes, probe2af){
  pAA_prior <- (1 - probe2af) ^2 # Prior genotype probability of AA
  pAB_prior <- 2 * probe2af * (1 - probe2af)
  pBB_prior <- probe2af ^2
  pAA_prior <- matrix(rep(pAA_prior, ncol(RAI)), ncol=ncol(RAI), dimnames=list(rownames(RAI), colnames(RAI)))
  pAB_prior <- matrix(rep(pAB_prior, ncol(RAI)), ncol=ncol(RAI), dimnames=list(rownames(RAI), colnames(RAI)))
  pBB_prior <- matrix(rep(pBB_prior, ncol(RAI)), ncol=ncol(RAI), dimnames=list(rownames(RAI), colnames(RAI)))
  shapes$mean <- shapes$shape1 / (shapes$shape1 + shapes$shape2)
  shapes <- as.matrix(shapes[order(shapes$mean),]) # row1 to row3: AA, AB, BB
  pD_AA <- apply(RAI, 1:2, function(x) dbeta(x, shapes[1, 1], shapes[1, 2])) # probability of data given genotype AA
  pD_AB <- apply(RAI, 1:2, function(x) dbeta(x, shapes[2, 1], shapes[2, 2]))
  pD_BB <- apply(RAI, 1:2, function(x) dbeta(x, shapes[3, 1], shapes[3, 2]))
  pD <- pD_AA * pAA_prior + pD_AB * pAB_prior + pD_BB * pBB_prior
  pAA <- pD_AA * pAA_prior / pD
  pAB <- pD_AB * pAB_prior / pD
  pBB <- pD_BB * pBB_prior / pD
  list(pAA=pAA, pAB=pAB, pBB=pBB)
}

#' Call genotypes based on EM algorithm and Bayesian approach
#' 
#' The Expectation–maximization (EM) algorithm is used to fit a mixture of three beta distributions representing the three genotypes (AA, AB, and BB). Then, the Bayesian approach is used to get the genotype probabilities, with AFs of the matched population (the 1000 Genomes Project, 1KGP) being used as priors.
#'
#' @param RAI A matrix of RAI (Ratio of Alternative allele Intensity) for probes. Provide probes as rows and samples as columns.
#' @param pop Population to be used to extract AFs. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @param maxiter Maximal number of iterations for the EM algorithm.
#' @param plotIter To plot the distribution of RAIs for each iteration.
#' @return  A list containing
#' \item{RAI}{Ratio of Alternative allele Intensity}
#' \item{fits}{Parameters of the mixture model}
#' \item{GP}{Posterior probabilities for the three genotypes}
#' @export
call_genotypes_bayesian <- function(RAI, pop, type, maxiter=50, plotIter=FALSE){
  
  #### deal with NAs
  
  # Initialize: classify all RAIs into three clusters
  RAI_plain <- RAI %>% as.data.frame %>%
    mutate(Probe=rownames(.)) %>%
    tidyr::gather(key="Sample", value="RAI", -Probe) %>% tibble
  RAI_plain$Cluster <- sapply(RAI_plain$RAI, function(x){if(x<.3){"Cluster1"}else if(x<.7){"Cluster2"}else{"Cluster3"}})
  
  # EM
  iterations <- list()
  iterations[[1]] <- iterate_em(RAI_plain)
  print(paste0("Iteration 1, -log(L) scaled by number of data points: ", paste(round(iterations[[1]]$fits$mLL, 3), collapse=" ")))
  i <- 2
  maxDiff <- Inf
  while(i<maxiter & maxDiff>0.1){
    iterations[[i]] <- iterate_em(iterations[[i-1]]$assignments)
    maxDiff = max(abs(iterations[[i]]$fits$mLL - iterations[[i-1]]$fits$mLL))
    print(paste0("Iteration ", i, ", -log(L) scaled by number of data points: ", paste(round(iterations[[i]]$fits$mLL, 3), collapse=" ")))
    i=i+1
  }
  finalClusters <- iterations[[i-1]]
  
  # Plots
  probe2af <- get_AF(RAI, pop, type)
  if(plotIter){plotIter_func(iterations, type)}

  # Get posterior genotype probabilities
  GP <- get_GP(RAI, finalClusters$fits[, c("shape1", "shape2")], probe2af)
  
  # return
  list(iterations=iterations, RAI=RAI, fits=finalClusters$fits, GP=GP)
}

#' Plot distribution of RAIs for each iteration
#'
#' @param iterations
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @export
plotIter_func <- function(iterations, type){
  iterations_assignments <- list.map(iterations, assignments)
  iterations_assignments <- do.call(
    rbind, 
    lapply(
      1:length(iterations_assignments), 
      function(i) tibble(iteration=paste0("Iteration ",i), iterations_assignments[[i]]))
  )
  iterations_fits <- list.map(iterations, fits)
  iterations_fits <- do.call(
    rbind, 
    lapply(
      1:length(iterations_fits), 
      function(i) tibble(iteration=paste0("Iteration ",i), iterations_fits[[i]]))
  ) 
  iterations_fits_density <- do.call(rbind, apply(
    iterations_fits, 1,
    function(x) tibble(
      iteration=x["iteration"], Cluster=x["Cluster"], index=0:100/100,
      density=as.numeric(x["prior"]) * dbeta(0:100/100, as.numeric(x["shape1"]), as.numeric(x["shape2"]))
    )))
  p <- ggplot(iterations_assignments) + 
    geom_histogram(aes(x=RAI, fill=Cluster), bins=100, alpha=.5) +
    geom_line(aes(x=index, y=density*30, color=Cluster), data=iterations_fits_density) +
    facet_wrap(~ iteration) + theme_bw() + labs(y="Counts or density")
  ggsave(filename=paste0("RAI_in_iteration.", type, ".pdf"), plot=p, width=12, height=12, units="in", scale=2)
}


