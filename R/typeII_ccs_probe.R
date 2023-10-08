
#' Call genotypes for Type II CCS probes
#' 
#' @param inData If input="raw", provide rgData here (Noob and dye-bias corrected signals produced by using `correct_noob_dye`). Otherwise, provide beta or M-value matrix here.
#' @param input Input data types. One of "raw", "beta", and "mval". If input is "beta" or "mval", please use probes as rows and samples as columns.
#' @param plotRAI If TRUE, plot distribution of RAIs.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param vcfName VCF file name. Only effective when vcf=TRUE.
#' @param a2 Call genotypes for specific alternative alleles. Use "AT" or "CG".
#' @param GP_cutoff When calculating missing rate, genotypes with the highest genotype probability < GP_cutoff will be treated as missing.
#' @param outlier_cutoff "max" or a number ranging from 0 to 1. If outlier_cutoff="max", genotypes with outlier probability larger than all of the three genotype probabilities will be set as missing. If outlier_cutoff is a number, genotypes with outlier probability > outlier_cutoff will be set as missing.
#' @param missing_cutoff Missing rate cutoff to filter variants. Note that for VCF output, variants with missing rate above the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with missing rate above the cutoff will be removed.
#' @param R2_cutoff_up,R2_cutoff_down R-square cutoffs to filter variants (Variants with R-square > R2_cutoff_up or < R2_cutoff_down should be removed). Note that for VCF output, variants with R-square outside this range will be marked in the `FILTER` column. For the returned dosage matrix, variants with R-square outside this range will be removed.
#' @param MAF_cutoff A MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with MAF below the cutoff will be removed.
#' @param HWE_cutoff HWE p value cutoff to filter variants. Note that for VCF output, variants with HWE p value below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with HWE p value below the cutoff will be removed.
#' @param train If TRUE, will fit the distribution of RAI (Ratio of Alternative allele Intensity) and filter probes by number of peaks. If FALSE, will use predefined probe list.
#' @param cpu Number of CPU. Only effective when train=TRUE.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL. Only probes with MAF of matching population > 0.01 will be kept. Only effective when train=TRUE.
#' @param maxiter Maximal number of iterations for the EM algorithm.
#' @param bayesian Use the Bayesian approach to calculate posterior genotype probabilities.
#' @param platform EPIC or 450K.
#' @param verbose Verbose mode: 0/1/2.
#' @return A list containing
#' \item{dosage}{A matrix of genotype calls. Variants with R2s, HWE p values, MAFs, or missing rates beyond the cutoffs are removed.}
#' \item{genotypes}{A list containing RAI, shapes of the mixed beta distributions, prior probabilities that the RAI values belong to one of the three genotypes, proportion of RAI values being outlier (U), and genotype probability (GP).}
#' \item{methyl_recalc}{Re-calculated methylation levels on reference alleles. A list containing shapes of the mixed beta distributions and true methylation level (pM) for each probe.}
#' @export
callGeno_typeII <- function(inData, input="raw", plotRAI=FALSE, vcf=FALSE, vcfName="genotypes.typeII_ccs_probe.vcf", a2="AT",
                            GP_cutoff=0.9, outlier_cutoff="max", missing_cutoff=0.1, 
                            R2_cutoff_up=1.1, R2_cutoff_down=0.75, MAF_cutoff=0.01, HWE_cutoff=1e-6, 
                            train=TRUE, cpu=1, pop="EAS", maxiter=50, bayesian=FALSE, platform="EPIC", verbose=1){
  if(!train & platform!="EPIC"){
    print("Error: train=FALSE only works with platform=EPIC.")
    return(NA)
  }
  if(!(pop %in% c("EAS", "AMR", "AFR", "EUR", "SAS", "ALL"))){
    stop("pop must be one of EAS, AMR, AFR, EUR, SAS, and ALL.")
  }
  tag_af <- paste0(pop, "_AF")
  if(platform=="EPIC"){
    if(a2=="AT"){
      data(probeInfo_typeII)
    }else if(a2=="CG"){
      data(probeInfo_typeII_CG); probeInfo_typeII <- probeInfo_typeII_CG
    }else{
      stop("Error: wrong a2 specified!")
    }
  }else{
    if(a2=="AT"){
      data(probeInfo_typeII_450K); probeInfo_typeII <- probeInfo_typeII_450K
    }else if(a2=="CG"){
      data(probeInfo_typeII_CG_450K); probeInfo_typeII <- probeInfo_typeII_CG_450K
    }else{
      stop("Error: wrong a2 specified!")
    }
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
  
  # get positions of the central modes
  if(train){
    mod <- getMod(beta, cpu=cpu)
  }else{
    mod <- probeInfo_typeII %>% dplyr::select(SNP, CpG, loc_pass, nmod, loc0, loc1, loc2)
    rownames(mod) <- mod$CpG
  }
  mod <- mod[mod$nmod > 1,]

  # calculate RAI
  beta <- beta[rownames(beta) %in% rownames(mod), ]
  if(a2=="AT"){
    mod$pM <- sapply(2 * mod$loc1, function(x) min(x, 1))
    RAI <- 1 - ( beta / matrix(rep(mod[rownames(beta), "pM", drop=TRUE], ncol(beta)), nrow=nrow(beta)) )
  }else if(a2=="CG"){
    mod$pM <- sapply(2 * mod$loc1 - 1, function(x) max(x, 0))
    pM_matrix <- matrix(rep(mod[rownames(beta), "pM", drop=TRUE], ncol(beta)), nrow=nrow(beta))
    RAI <- ( beta - pM_matrix ) / ( 1 - pM_matrix )
  }else{
    stop("Error: wrong a2 specified!")
  }
  RAI[RAI < 0.01] <- 0.01 # if set to 0 or 1, it will fail in fitting beta distribution as GP will be NA in call_genotypes.R:: GP <<- GP / tmp
  RAI[RAI > 0.99] <- 0.99
  
  # filter probes based on peak density and positions.
  if(train){
    mod <- getMod(RAI, cpu=cpu)
    RAI <- RAI[dplyr::filter(mod, loc_pass==TRUE, nmod > 1)$CpG,]
  }else{
    mod <- dplyr::filter(probeInfo_typeII, loc_pass==TRUE, nmod > 1) %>% 
      dplyr::select(SNP, CpG, loc_pass, nmod, loc0, loc1, loc2)
    rownames(mod) <- mod$CpG
    RAI <- RAI[mod$CpG,]
  }
  
  # call genotypes
  genotypes <- call_genotypes(RAI, pop=pop, type="typeII_ccs_probe", maxiter=maxiter, 
                                       bayesian=bayesian, platform=platform, verbose=verbose)
  if(plotRAI){plot_RAI_distribution(genotypes, type="typeII_ccs_probe")}
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, a2=a2,
                             GP_cutoff=GP_cutoff, outlier_cutoff=outlier_cutoff, missing_cutoff=missing_cutoff,
                             R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, 
                             MAF_cutoff=MAF_cutoff, HWE_cutoff=HWE_cutoff, 
                             type="typeII_ccs_probe", pop=pop, plotAF=FALSE, platform=platform)
  list(dosage=dosage, mod=mod, genotypes=genotypes)
}

