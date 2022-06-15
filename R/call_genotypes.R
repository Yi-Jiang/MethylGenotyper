#' Call genotypes (adapted from `ewastools`)
#' 
#' A mixture model (3 Beta distributions, 1 uniform distribution for outliers) is fitted to probes. The model parameters are learnt via Expectation-Maximization (EM) algorithm. We used a "two-step" strategy to call genotypes, which is run by fitting the mixture model twice. For the first step, we fit the mixture model, call genotypes, and filter out the probes deviating Hardy-Weinberg Equilibrium (HWE, P < 1E-6). For the second step, we fit the model and call genotypes again, using the pre-filtered data.
#'
#' @param RAI A matrix of RAI (Ratio of Alternative allele Intensity) for probes. Provide probes as rows and samples as columns.
#' @param learn If learn=TRUE, will learn dataset-specific parameters using the EM algorithm. The genotype calling process will be iterated two times, with the first iteration for selecting probes and the second iteration for re-training.
#' @param maxiter Maximal number of iterations of the EM algorithm learning the mixture model.
#' @return  A list containing
#' \item{par}{Parameters of the mixture model}
#' \item{loglik}{Log-likelihood in each iteration of the EM algorithm}
#' \item{outliers}{A-posteriori probability of SNP being an outlier}
#' \item{gamma}{A-posteriori probabilities for each of the three genotypes}
#' @export
callGeno <- function(RAI, learn=FALSE, maxiter=50){
  probes = rownames(RAI)
  probe_marks = rep(probes, times=ncol(RAI)) # Mark probe positions
  RAI_vec = RAI
  dim(RAI_vec) = NULL
  
  # Drop NAs to be able to compute likelihoods, but keep score of which entries to reinsert them again later
  NAs = is.na(RAI_vec)
  RAI_vec_hq = RAI_vec[!NAs]
  n = length(RAI_vec_hq)

  if(learn==FALSE){ # Use predefined model parameters (might work better if training set is small or contains many outliers)
    # Class probability for outliers
    alpha = 0.06646095
    
    # Class probabilities for homozygous and heterozygous genotypes
    pi = c(0.2818387, 0.4330363, 0.2851250)
    
    # Beta distribution parameters
    shapes1 = c(2.206479, 80.830012, 40.640821)
    shapes2 = c(38.043029, 84.411900, 3.315509)
    
    # Uniform distribution representing outliers
    p = 1
    
    loglik = NULL
    gamma = cbind(
      pi[1] * dbeta(RAI_vec_hq, shape1=shapes1[1], shape2=shapes2[1]),
      pi[2] * dbeta(RAI_vec_hq, shape1=shapes1[2], shape2=shapes2[2]),
      pi[3] * dbeta(RAI_vec_hq, shape1=shapes1[3], shape2=shapes2[3])
    )
    gamma = (1-alpha) * gamma
    tmp = rowSums(gamma)
    gamma = gamma/tmp
    outliers = (alpha*p) / ((alpha*p) + tmp)
  }else{ # Learn dataset-specific model parameters using the EM algorithm
    # Initialize
    alpha = 1e-2 # Class probability for outliers
    outliers = rep(alpha,times=n)
    pi = c(1/3, 1/3, 1/3) # Class probabilities for homozygous and heterozygous genotypes
    shapes1 = c(10, 80, 80) # Beta distribution parameters
    shapes2 = c(80, 80, 10)
    p = 1 # Uniform distribution representing outliers
    gamma = NA
    
    e_step = function(){
      gamma = cbind(
        pi[1] * dbeta(RAI_vec_hq, shape1=shapes1[1], shape2=shapes2[1]),
        pi[2] * dbeta(RAI_vec_hq, shape1=shapes1[2], shape2=shapes2[2]),
        pi[3] * dbeta(RAI_vec_hq, shape1=shapes1[3], shape2=shapes2[3])
      )
      gamma = (1-alpha) * gamma # gamma is a (m+n)*3 matrix, where m is number of probes, n is number of samples, columns are three genotypes.
      tmp = rowSums(gamma)
      gamma <<- gamma/tmp
      outliers <<- (alpha*p) / (alpha*p + tmp)
      loglik = sum(log(alpha*p + tmp))
      loglik
    }
    
    m_step = function(){
      gamma = gamma * (1-outliers)
      # MLE
      s1 = eBeta(RAI_vec_hq, gamma[,1])
      s2 = eBeta(RAI_vec_hq, gamma[,2])
      s3 = eBeta(RAI_vec_hq, gamma[,3])
      shapes1 <<- c(s1$shape1, s2$shape1, s3$shape1)
      shapes2 <<- c(s1$shape2, s2$shape2, s3$shape2)
      # MLE of class priors
      pi = apply(gamma, 2, sum)
      pi <<- pi/sum(pi)
      alpha <<- sum(outliers)/n
      invisible(NULL)
    }
    
    runEM = function(){
      loglik = rep(NA_real_, times=maxiter)
      loglik[1] = e_step()
      i = 2; gain=Inf;
      while(i<maxiter & gain>1e-4){ # stop if maxiter reached or improvement is below threshold
        m_step()
        loglik[i] = e_step()
        gain = loglik[i]-loglik[i-1]
        i=i+1
      }
      loglik=loglik[1:(i-1)]
      loglik
    }
    
    # First iteration
    loglik <- runEM()
    
    # Remove probes of low quality
    fail_1st <- rep(FALSE, length(probe_marks))
    fail_1st_probes <- c()
    for(probe in probes){
      gamma0 <- gamma[probe_marks[!NAs]==probe,] # gamma values for a specific probe. row: sample, col: three genotype probabilities.
      colnames(gamma0) <- c("0/0", "0/1", "1/1")
      hardgeno <- apply(gamma0, 1, function(x) names(which.max(x)))
      hardgeno_sum <- matrix(
        c(sum(hardgeno=="0/0"), sum(hardgeno=="0/1"), sum(hardgeno=="1/1")),
        nrow=1, ncol=3, dimnames=list(NULL, c("MM", "MN", "NN"))
      )
      hwe <- suppressWarnings(HWChisqMat(hardgeno_sum, cc=0, verbose=FALSE))
      #dosage <- gamma0[,2] + 2 * gamma0[,3]
      #AF <- mean(dosage) / 2
      #R2 <- var(dosage) / (2 * AF * (1 - AF))
      #if(is.na(R2) | R2 < 0.5 | R2 > 1.3){
      if(hwe$pvalvec < 1e-6){
        fail_1st[probe_marks==probe] <- TRUE
        fail_1st_probes <- c(fail_1st_probes, probe)
      }
    }
    RAI_vec_hq = RAI_vec[!NAs & !fail_1st]
    n = length(RAI_vec_hq)
    
    # Second iteration
    loglik <- runEM()
  }
  
  ## Re-insert missing values
  tmp = rep(NA_real_, times=length(RAI))
  tmp[!NAs & !fail_1st] = outliers
  outliers = tmp
  dim(outliers) = dim(RAI)
  gamma = lapply(1:3,function(k){
    tmp = rep(NA_real_, times=length(RAI))
    tmp[!NAs & !fail_1st] = gamma[,k]
    dim(tmp) = dim(RAI)
    tmp
  })
  
  return(list(
    snps=RAI,
    outliers=outliers, # m*n, where m is number of probes, n is number of samples.
    gamma=gamma, # list of 3, indicating genotype probabilities. Each element is a m*n matrix.
    par=list(pi=pi, shapes1=shapes1, shapes2=shapes2, alpha=alpha),
    loglik=loglik,
    fail_1st_probes=fail_1st_probes # probes failed the first iteration.
  ))
}

#' eBeta (adapted from `ewastools`)
#' @export
eBeta = function(x,w){
  
  # Beta distribution parameter estimation
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