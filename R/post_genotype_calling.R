
#' Get sample relationships
#' 
#' @param Phi A vector of Phi.
#' @return A vector of sample relationships.
#' @export
getRelation <- function(phi) {
  a <- ifelse(phi>=0.5^1.5, "MZ",
    ifelse(phi>=0.5^2.5, "1st",
      ifelse(phi>=0.5^3.5, "2nd",
        ifelse(phi>=0.5^4.5, "3rd", "UN")
      )
    )
  )
  a <- factor(a, levels=c("MZ", "1st", "2nd", "3rd", "UN"))
  a
}

#' Get kinship coefficients using the SEEKIN estimator
#' 
#' For homogeneous samples. Only SNPs with missing rate < 10% were used.
#' 
#' @param dosage Dosage genotypes.
#' @return A data frame containing kinship coefficient (Phi) and sample relationships between each two samples.
#' @export
getKinship <- function(dosage){
  missing <- apply(dosage, 1, function(x) sum(is.na(x))) / ncol(dosage)
  nValue <- apply(dosage, 1, function(x) length(unique(x[!is.na(x)])))
  dosage <- dosage[missing < 0.1 & nValue >1,]
  AF <- rowMeans(dosage, na.rm=TRUE) / 2
  R2 <- apply(dosage, 1, function(x) var(x, na.rm=TRUE)) / (2 * AF * (1 - AF))
  R2[R2>1] <- 1
  norm <- sweep(dosage, 1, 2*AF, FUN="-")
  norm[is.na(norm)] <- 0
  GRM <- t(norm) %*% norm / sum(2 * AF * (1 - AF) * R2^2)
  GRM[!(lower.tri(GRM))] <- NA
  GRM <- as.data.frame(GRM) %>% mutate(IID1=rownames(.)) %>%
    tidyr::gather(key="IID2", value="Relatedness", -IID1) %>%
    dplyr::filter(!is.na(Relatedness)) %>%
    mutate(Phi=Relatedness/2, Relation=getRelation(Phi))
  GRM
}

#' Get kinship coefficients using the SEEKIN-het estimator
#' 
#' For heterogeneous samples. Only SNPs with missing rate < 10% were used.
#' 
#' @param dosage A matrix of genotype calls. Provide probes as rows and samples as columns.
#' @param indAF A matrix of individual-specific AFs. Provide probes as rows and samples as columns.
#' @return A data frame containing kinship coefficient (Phi) and sample relationships between each two samples.
#' @export
getKinship_het <- function(dosage, indAF){
  probes <- intersect(rownames(dosage), rownames(indAF))
  samples <- intersect(colnames(dosage), colnames(indAF))
  dosage <- dosage[probes, samples]
  indAF <- indAF[probes, samples]
  missing <- apply(dosage, 1, function(x) sum(is.na(x))) / ncol(dosage)
  nValue <- apply(dosage, 1, function(x) length(unique(x[!is.na(x)])))
  dosage <- dosage[missing < 0.1 & nValue >1,]
  AF <- rowMeans(dosage, na.rm=TRUE)/2
  R2 <- apply(dosage, 1, function(x) var(x, na.rm=TRUE))/(2 * AF * (1 - AF))
  R2[R2 > 1] <- 1
  mu <- AF+R2*(indAF-rowMeans(indAF))
  norm <- dosage-2*mu
  norm[is.na(norm)] <- 0
  phat <- sqrt(2*indAF*(1-indAF)*R2^2)
  GRM <- t(norm)%*%norm/(t(phat) %*% phat)
  GRM[!(lower.tri(GRM))] <- NA
  GRM <- as.data.frame(GRM) %>% mutate(IID1=rownames(.)) %>%
    tidyr::gather(key="IID2", value="Relatedness", -IID1) %>%
    dplyr::filter(!is.na(Relatedness)) %>%
    mutate(Phi=Relatedness/2, Relation=getRelation(Phi))
  GRM
}

#' Predict sample contamination according to inbreeding coefficients
#' 
#' Samples with inbreeding coefficients beyond 3 SDs of the mean are considered as contaminated.
#' 
#' @param F A vector of inbreeding coefficients.
#' @return A vector of whether samples are contaminated (outlier) or not contaminated (normal).
#' @export
predContamination <- function(x) {
  upper <- mean(x)+3*sd(x)
  lower <- mean(x)-3*sd(x)
  ifelse(x>upper, "outlier",
         ifelse(x<lower, "outlier", "normal")
  )
}

#' Calculate inbreeding coefficients according to Genotypes
#' 
#' @param genotype Genotype matrix, with each row indicates a SNP and each column indicates a sample. Dosage genotypes are forced to binary (0, 1, or 2) by using the `round` function.
#' @return A data frame of inbreeding coefficients and sample contamination status.
#' @export
getInbreed <- function(dosage){
  binGeno <- round(dosage, 0) # force to 0/1/2
  pHet_exp <- apply(
    binGeno, 1, 
    function(x){
      m <- mean(x, na.rm=TRUE)
      m*(1-m/2)
    }
  )
  samples <- colnames(binGeno)
  nHet_obs <- c()
  nHet_exp <- c()
  for (sp in samples){
    idx <- which(!is.na(binGeno[,sp]))
    nHet_obs[sp] <- sum(binGeno[idx, sp]==1)
    nHet_exp[sp] <- sum(pHet_exp[idx])
  }
  inbreed <- tibble(IID = samples, 
                    nHet_obs = nHet_obs[samples],
                    nHet_exp = nHet_exp[samples],
                    F_pred = 1 - nHet_obs[samples] / nHet_exp[samples],
                    contamination = predContamination(1 - nHet_obs[samples] / nHet_exp[samples]))
  inbreed
}
