
#' PCA and Procrustes analysis
#'
#' @param studyGeno A matrix of genotypes of study samples. Provide probes as rows and samples as columns. Include all SNP probes, type I CCS probes, and type II CCS probes if available.
#' @param plotPCA To plot the projection of study samples in reference ancestry space.
#' @param cpu Number of CPU.
#' @return A list containing
#' \item{refPC}{Top PCs in the reference}
#' \item{studyPC}{Top PCs in study samples}
#' @export
projection <- function(studyGeno, plotPCA=TRUE, cpu=1){
  ## Filter SNPs for PCA: MAF>0.1 and R2>0.9
  data(cpg2snp)
  AF <- rowMeans(studyGeno) / 2
  R2 <- apply(studyGeno, 1, var) / (2 * AF * (1 - AF))
  studyGeno <- studyGeno[AF>0.1 & R2>0.9,]
  rownames(studyGeno) <- cpg2snp[rownames(studyGeno)]
  
  ## PCA and Procrustes analysis
  data(refGeno_1KGP3)
  pc <- trace(refGeno_1KGP3, studyGeno, cpu=cpu)
  
  ## Plot PCA
  if(plotPCA){
    plotPCA(pc$refPC, pc$studyPC)
  }
  
  pc
}

#' Recalibrate genotypes for samples of mixed population
#'
#' @param genotypes A list returned by either `callGeno_snp`, `callGeno_typeI`, or `callGeno_typeII` function.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @param refPC Top PCs in the reference
#' @param studyPC Top PCs in study samples
#' @return A list of recalibrated genotypes containing
#' \item{dosage}{A matrix of genotype calls}
#' \item{genotypes}{A list containing RAI, fits, and Genotype probabilities}
#' @export
recal_Geno <- function(genotypes, type, refPC, studyPC){
  data(cpg2snp)
  data(snp2cpg)
  ## Model genotypes of the reference individuals as a linear function of PCs
  betas <- list()
  for(snp in cpg2snp[rownames(genotypes$genotypes$RAI)]){
    betas[[snp]] <- coefficients(lm(refGeno_1KGP3[snp,] ~ refPC))
  }
  betas <- do.call(rbind, betas)
  colnames(betas) <- c("Intercept", paste0("RefPC", 1:4))
  
  ## Calculate individual-specific AFs
  studyPC <- cbind(Intercept=1, studyPC)
  indAF <- betas %*% t(studyPC) / 2
  indAF[indAF<0.001] <- 0.001 # constrain AFs to avoid out of boundary values
  indAF[indAF>0.999] <- 0.999
  rownames(indAF) <- snp2cpg[rownames(indAF)]
  
  ## Recalibrate posterior genotype probabilities
  GP <- get_GP(genotypes$genotypes$RAI, genotypes$genotypes$fits[, c("shape1", "shape2")], indAF)
  genotypes_recal <- list(genotypes = genotypes$genotypes)
  genotypes_recal$genotypes$GP <- GP
  genotypes_recal$genotypes$RAI <- genotypes_recal$genotypes$RAI[rownames(GP$pAA), colnames(GP$pAA)]
  genotypes_recal$dosage <- format_genotypes(genotypes_recal$genotypes, vcf=T, 
                                             vcfName=paste0("genotypes.recal.", type, ".vcf"),
                                             R2_cutoff_up=1.1, R2_cutoff_down=0.75, MAF_cutoff=0.01,
                                             type=type, plotAF=FALSE)
  genotypes_recal
}

#' To plot the projection of study samples in reference ancestry space
#'
#' @param refPC Top PCs in the reference
#' @param studyPC Top PCs in study samples
#' @export
plotPCA <- function(refPC, studyPC){
  data(sam2pop)
  refPC <- data.frame(refPC, popID=sam2pop[rownames(refPC)])
  studyPC <- data.frame(studyPC, popID="Study")
  p1 <- ggplot(refPC) +
    geom_point(aes(x=PC1, y=PC2, color=popID), size=3, alpha=0.6) +
    scale_color_brewer(palette="Set2")+
    geom_point(aes(x=PC1, y=PC2, shape=popID), size=3, alpha=0.6, data=studyPC) +
    scale_shape_manual(values=c("Study"=1))+
    theme_bw()
  p2 <- ggplot(refPC) +
    geom_point(aes(x=PC3, y=PC4, color=popID), size=3, alpha=0.6) +
    scale_color_brewer(palette="Set2")+
    geom_point(aes(x=PC3, y=PC4, shape=popID), size=3, alpha=0.6, data=studyPC) +
    scale_shape_manual(values=c("Study"=1))+
    theme_bw()
  p <- ggarrange(p1, p2, nrow=1, ncol=2)
  ggsave(filename="trace.pca.pdf", plot=p, width=7, height=3, units="in", scale=2)
}
