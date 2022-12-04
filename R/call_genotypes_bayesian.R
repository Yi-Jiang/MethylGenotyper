
#' The Expectation–maximization (EM) algorithm: M step
#' 
#' @param assignments A data frame containing:
#' \item{Probe}{Probe ID}
#' \item{Sample}{Sample ID}
#' \item{RAI}{Ratio of Alternative allele Intensity}
#' \item{Cluster}{Cluster IDs}
#' @return A list of two elements:
#' \item{parameters}{A data frame of: 1) cluster ID, 2) shape1, 3) shape2, 4) probability of being an outlier, and 5) prior probabilitiy of the genotype being AA,AB,or BB}
#' \item{mLL_norm}{Minus log-likelihood scaled by number of data points}
#' @export
m_step <- function(x){
  N <- nrow(x)
  priors <- table(x$Cluster) / N
  minuslogL <- function(shape1_1, shape2_1, shape1_2, shape2_2, shape1_3, shape2_3, U){
    #write.table(paste(c(shape1_1, shape2_1, shape1_2, shape2_2, shape1_3, shape2_3, U), collapse="\t"), file="tmp", append=T, quote=F, row.names=F, col.names=F)
    -sum(
      log(U + (1-U) * rowSums(cbind(
        priors["Cluster1"] * dbeta(x$RAI, shape1_1, shape2_1, log=F),
        priors["Cluster2"] * dbeta(x$RAI, shape1_2, shape2_2, log=F),
        priors["Cluster3"] * dbeta(x$RAI, shape1_3, shape2_3, log=F)
      )))
    )
  }
  m <- stats4::mle(
    minuslogL, 
    start=list(
      shape1_1=5,  shape2_1=60,
      shape1_2=30, shape2_2=30,
      shape1_3=60, shape2_3=5, U=0.01), 
    method="L-BFGS-B", 
    lower = c(rep(0.1, 6), 0.001),
    upper = c(rep(150, 6), 0.3), 
    control = list(ndeps = c(rep(0.1, 6), 0.001), trace=0)
  )
  s <- stats4::coef(m)
  parameters <- tibble(
    Cluster = paste0("Cluster", 1:3),
    shape1 = c(s["shape1_1"], s["shape1_2"], s["shape1_3"]),
    shape2 = c(s["shape2_1"], s["shape2_2"], s["shape2_3"]),
    U = s["U"],
    prior = priors
  )
  list(parameters = parameters, mLL_norm = m@min/N)
}

#' The Expectation–maximization (EM) algorithm: E step
#' 
#' @param assignments A data frame containing:
#' \item{Probe}{Probe ID}
#' \item{Sample}{Sample ID}
#' \item{RAI}{Ratio of Alternative allele Intensity}
#' \item{Cluster}{Cluster IDs}
#' @param fits A list produced in m_step
#' @return A data frame with re-assigned cluster IDs
#' @export
e_step <- function(assignments, fits){
  assignments <- assignments %>%
    dplyr::select(Probe:RAI) %>%
    crossing(fits$parameters) %>%
    mutate(likelihood = prior * dbeta(RAI, shape1, shape2)) %>%
    group_by(Probe, Sample) %>%
    top_n(1, likelihood) %>%
    ungroup()
  assignments
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
#' @param AF A MxN matrix of AFs. Provide SNPs as rows and samples as columns.
#' @return  A list containing
#' \item{pAA}{Posterior genotype probability of AA}
#' \item{pAB}{Posterior genotype probability of AB}
#' \item{pBB}{Posterior genotype probability of BB}
#' @export
get_GP <- function(RAI, shapes, AF){
  probes <- intersect(rownames(RAI), rownames(AF))
  samples <- intersect(colnames(RAI), colnames(AF))
  RAI <- RAI[probes, samples]
  AF <- AF[probes, samples]
  pAA_prior <- (1 - AF) ^2 # Prior genotype probability of AA
  pAB_prior <- 2 * AF * (1 - AF)
  pBB_prior <- AF ^2
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
#' The Expectation–maximization (EM) algorithm is used to fit a mixture of three beta distributions representing the three genotypes (AA, AB, and BB). Then, the Bayesian approach is used to get the genotype probabilities, with AFs of the matched population (the 1000 Genomes Project, 1KGP) being used to infer priors.
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
  # Initialize: classify all RAIs into three clusters
  assignments <- RAI %>% as.data.frame %>%
    mutate(Probe=rownames(.)) %>%
    tidyr::gather(key="Sample", value="RAI", -Probe) %>% tibble
  #assignments <- assignments[sample(1:nrow(assignments), 100000),]
  assignments$Cluster <- sapply(assignments$RAI, function(x){if(x<.3){"Cluster1"}else if(x<.7){"Cluster2"}else{"Cluster3"}})
  
  # EM
  iterations <- list()
  maxDiff <- Inf
  i <- 1
  while(i < maxiter & maxDiff > 1e-4){
    fits <- m_step(assignments)
    assignments <- e_step(assignments, fits)
    iterations[[i]] <- fits
    #iterations[[i]]$assignments <- assignments
    print(paste0("Iteration ", i, ", -log(L) scaled by number of data points: ", round(fits$mLL_norm, 6)))
    if(i>1){
      maxDiff = abs(iterations[[i]]$mLL_norm - iterations[[i-1]]$mLL_norm)
    }
    i=i+1
  }
  finalClusters <- iterations[[i-1]]
  
  # Plots
  if(plotIter){plotIter_func(iterations, type)}

  # Get posterior genotype probabilities
  probe2af <- get_AF(pop, type)
  AF <- matrix(rep(probe2af[rownames(RAI)], ncol(RAI)), ncol=ncol(RAI), dimnames=list(rownames(RAI), colnames(RAI)))
  GP <- get_GP(RAI, finalClusters$fits[, c("shape1", "shape2")], AF)
  
  # return
  list(RAI=RAI, fits=finalClusters$fits, iterations=iterations, GP=GP)
}

#' Plot distribution of RAIs for each iteration
#'
#' @param iterations EM iteration results
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


