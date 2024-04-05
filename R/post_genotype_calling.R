
#' Get sample relationships
#' 
#' @param phi A vector of kinship coefficient (Phi).
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

#' Get kinship coefficients and inbreeding coefficients using the SEEKIN estimator
#' 
#' Only SNPs with missing rate < 10% were used.
#' 
#' @param dosage A matrix of genotype calls. Provide probes as rows and samples as columns.
#' @return A list containing
#' \item{kinship}{A data frame containing kinship coefficient (Phi) and sample relationships between each two samples.}
#' \item{inbreed}{A vector of inbreeding coefficients.}
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
  inbreed <- GRM[row(GRM)==col(GRM)] - 1 # inbreeding coefficients = 2 * Phi - 1
  names(inbreed) <- colnames(dosage)
  GRM[!(lower.tri(GRM))] <- NA
  kinship <- as.data.frame(GRM) %>% mutate(IID1=rownames(.)) %>%
    tidyr::gather(key="IID2", value="Relatedness", -IID1) %>%
    dplyr::filter(!is.na(Relatedness)) %>%
    mutate(Phi=Relatedness/2, Relation=getRelation(Phi))
  list(kinship=kinship, inbreed=inbreed)
}

#' Get kinship coefficients using the SEEKIN-het estimator
#' 
#' Only SNPs with missing rate < 10% were used.
#' 
#' @param dosage A matrix of genotype calls. Provide probes as rows and samples as columns.
#' @param indAF A matrix of individual-specific AFs. Provide probes as rows and samples as columns.
#' @return A data frame containing kinship coefficient (Phi) and sample relationships between each two samples.
#' @export
getKinship_het <- function(dosage, indAF){
  probes <- intersect(rownames(dosage), rownames(indAF))
  samples <- intersect(colnames(dosage), colnames(indAF))
  dosage <- dosage[probes, samples]
  missing <- apply(dosage, 1, function(x) sum(is.na(x))) / ncol(dosage)
  nValue <- apply(dosage, 1, function(x) length(unique(x[!is.na(x)])))
  dosage <- dosage[missing < 0.1 & nValue >1,]
  indAF <- indAF[rownames(dosage), samples]
  AF <- rowMeans(dosage, na.rm=TRUE)/2
  R2 <- apply(dosage, 1, function(x) var(x, na.rm=TRUE))/(2 * AF * (1 - AF))
  R2[R2 > 1] <- 1
  mu <- AF+R2*(indAF-rowMeans(indAF))
  norm <- dosage-2*mu
  norm[is.na(norm)] <- 0
  phat <- sqrt(2*indAF*(1-indAF)*R2^2)
  GRM <- t(norm)%*%norm/(t(phat) %*% phat)
  GRM[!(lower.tri(GRM))] <- NA
  kinship <- as.data.frame(GRM) %>% mutate(IID1=rownames(.)) %>%
    tidyr::gather(key="IID2", value="Relatedness", -IID1) %>%
    dplyr::filter(!is.na(Relatedness)) %>%
    mutate(Phi=Relatedness/2, Relation=getRelation(Phi))
  kinship
}
