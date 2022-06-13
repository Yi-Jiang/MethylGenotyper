
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

#' Get kinship coefficients and sample relationships
#' 
#' @param dosage Dosage genotypes.
#' @return A data frame containing kinship coefficient (Phi) and sample relationships between each two samples.
#' @export
getKinship <- function(dosage){
  AF <- rowMeans(dosage) / 2
  R2 <- apply(dosage, 1, var) / (2 * AF * (1 - AF))
  R2[R2>1] <- 1
  norm <- sweep(dosage, 1, 2*AF, FUN="-")
  GRM <- t(norm) %*% norm / sum(2 * AF * (1 - AF) * R2^2)
  GRM[!(lower.tri(GRM))] <- NA
  GRM <- as.data.frame(GRM) 
  GRM$IID1 <- rownames(GRM)
  GRM <- as.data.frame(melt(GRM))
  colnames(GRM) <- c("IID1", "IID2", "Relatedness")
  GRM <- GRM[!is.na(GRM$Relatedness),]
  GRM$Phi <- GRM$Relatedness / 2
  GRM$Relation <- getRelation(GRM$Phi)
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
  nHet_obs <- apply(binGeno, 2, function(x) sum(x==1))
  nHet_exp <- sum(apply(binGeno, 1, function(x) mean(x)*(1-mean(x)/2)))
  inbreed <- tibble(IID = colnames(binGeno), 
                    nHet_obs = nHet_obs,
                    nHet_exp = nHet_exp,
                    F_pred = 1 - nHet_obs / nHet_exp,
                    contamination = predContamination(1 - nHet_obs / nHet_exp))
  inbreed
}
