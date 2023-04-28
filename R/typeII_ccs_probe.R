
#' Call genotypes for Type II CCS probes
#' 
#' @param inData If input="raw", provide rgData here (Noob and dye-bias corrected signals produced by using `correct_noob_dye`). Otherwise, provide beta or M-value matrix here.
#' @param input Input data types. One of "raw", "beta", and "mval". If input is "beta" or "mval", please use probes as rows and samples as columns.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param vcfName VCF file name. Only effective when vcf=TRUE.
#' @param R2_cutoff_up,R2_cutoff_down R-square cutoffs to filter variants (Variants with R-square > R2_cutoff_up or < R2_cutoff_down should be removed). Note that for VCF output, variants with R-square outside this range will be marked in the `FILTER` column. For the returned matrix, variants with R-square outside this range will be removed.
#' @param MAF_cutoff A MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned matrix, variants with MAF below the cutoff will be removed.
#' @param train If TRUE, will fit the distribution of RAI (Ratio of Alternative allele Intensity) and filter probes by number of peaks. If FALSE, will use predefined probe list.
#' @param cpu Number of CPU. Only effective when train=TRUE.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL. Only probes with MAF of matching population > 0.01 will be kept. Only effective when train=TRUE.
#' @param maxiter Maximal number of iterations for the EM algorithm.
#' @param bayesian Use the Bayesian approach to calculate posterior genotype probabilities.
#' @param platform EPIC or 450K.
#' @param verbose Verbose mode: 0/1/2.
#' @return A list containing
#' \item{dosage}{A matrix of genotype calls. Variants with R2 or MAF beyond the cutoffs are removed. Genotypes with genotype quality (GQ) < 20 will be marked as NA.}
#' \item{genotypes}{A list containing RAI, shapes of the mixed beta distributions, prior probabilities that the RAI values belong to one of the three genotypes, proportion of RAI values being outlier (U), genotype probability (GP), Phred-scaled genotype likelihood (PL), and genotype quality (GQ).}
#' \item{methyl_recalc}{Re-calculated methylation levels on reference alleles. A list containing shapes of the mixed beta distributions and true methylation level (pM) for each probe.}
#' @export
callGeno_typeII <- function(inData, input="raw", plotBeta=FALSE, vcf=FALSE, vcfName="genotypes.typeII_ccs_probe.vcf", 
                            R2_cutoff_up=1.1, R2_cutoff_down=0.75, MAF_cutoff=0.01, train=TRUE, cpu=1, pop="EAS", 
                            maxiter=50, bayesian=TRUE, platform="EPIC", verbose=1){
  if(!train & platform!="EPIC"){
    print("Error: train=FALSE only works with platform=EPIC.")
    return(NA)
  }
  if(!(pop %in% c("EAS", "AMR", "AFR", "EUR", "SAS", "ALL"))){
    stop("pop must be one of EAS, AMR, AFR, EUR, SAS, and ALL.")
  }
  tag_af <- paste0(pop, "_AF")
  if(platform=="EPIC"){
    data(probeInfo_typeII)
  }else{
    data(probeInfo_typeII_450K); probeInfo_typeII <- probeInfo_typeII_450K
  }
  
  # calculate beta values
  if(input=="raw"){
    cg <- rownames(inData[["AR"]])
    cg_IIR <- cg[cg %in% probeInfo_typeII[probeInfo_typeII[,"Group"]=="IIR" & probeInfo_typeII[,tag_af]>0.01 & probeInfo_typeII[,tag_af]<0.99, "CpG"]] # Type II, Alt allele match Red
    beta <- inData[["AG"]][cg_IIR,] / (inData[["AG"]][cg_IIR,] + inData[["AR"]][cg_IIR,])
  }else if(input=="beta"){
    beta <- inData[rownames(inData) %in% probeInfo_typeII[probeInfo_typeII[,tag_af]>0.01 & probeInfo_typeII[,tag_af]<0.99, "CpG"], ]
  }else if(input=="mval"){
    beta <- mval2beta(inData[rownames(inData) %in% probeInfo_typeII[probeInfo_typeII[,tag_af]>0.01 & probeInfo_typeII[,tag_af]<0.99, "CpG"], ])
  }else{
    print("Error: Input data type must be one of raw, beta, and mval.")
    return(NA)
  }
  
  # calculate RAI
  print(paste(Sys.time(), "Calculating RAI values for type II CCS probes."))
  nCpG <- nrow(beta)
  pM <- c()
  shapes <- list()
  for(i in 1:nrow(beta)){
    CpG <- rownames(beta)[i]
    finalClusters <- fit_beta_em(beta[i,,drop=F], maxiter=maxiter, verbose=verbose-2)
    shapes[[CpG]] <- c(cluster1_shape1=finalClusters$shapes["Cluster1", "shape1"], cluster1_shape2=finalClusters$shapes["Cluster1", "shape2"],
                       cluster2_shape1=finalClusters$shapes["Cluster2", "shape1"], cluster2_shape2=finalClusters$shapes["Cluster2", "shape2"],
                       cluster3_shape1=finalClusters$shapes["Cluster3", "shape1"], cluster3_shape2=finalClusters$shapes["Cluster3", "shape2"],
                       U=finalClusters$U)
    cluster2 <- finalClusters$shapes["Cluster2", c("shape1", "shape2"), drop=T]
    cluster2_mean <- cluster2$shape1 / (cluster2$shape1 + cluster2$shape2)
    if(is.na(cluster2_mean)){cluster2_mean <- 0.5} # it's NA when the shapes are Inf.
    pM[CpG] <- min(2 * cluster2_mean, 1)
    if(verbose>=2){print(paste("Estimating true methylation level for Type II CCS probe", CpG, paste0("(", i, "/", nCpG, "):"), "pM =", pM[CpG]))}
  }
  RAI <- 1 - ( beta / matrix(rep(pM, ncol(beta)), nrow=nCpG) )
  RAI[RAI < 0.01] <- 0.01 # if set to zero, it will fail in fitting beta distribution as GP will be NA in call_genotypes_bayesian.R:: GP <<- GP / tmp
  RAI[RAI > 1] <- 1
  methyl_recalc <- list(shapes = do.call(rbind, shapes), trueMethyl = pM)

  # filter probes based on peak density and positions.
  if(train){
    mod <- getMod(RAI, cpu=cpu)
    RAI <- RAI[dplyr::filter(mod, h_0.1==TRUE, loc_pass==TRUE)$Name,]
  }else{
    RAI <- RAI[rownames(RAI) %in% dplyr::filter(probeInfo_typeII, h_0.1==TRUE, loc_pass==TRUE)$CpG,]
  }
  
  # call genotypes
  genotypes <- call_genotypes_bayesian(RAI, pop=pop, type="typeII_ccs_probe", maxiter=maxiter, 
                                       bayesian=bayesian, platform=platform, verbose=verbose)
  if(plotBeta){plot_beta_distribution(genotypes, type="typeII_ccs_probe")}
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, 
                             R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, 
                             MAF_cutoff=MAF_cutoff, type="typeII_ccs_probe", pop=pop, plotAF=FALSE, platform=platform)
  list(dosage=dosage, methyl_recalc=methyl_recalc, genotypes=genotypes)
}

