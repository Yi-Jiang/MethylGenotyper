
#' Call genotypes for Type II CCS probes
#' 
#' @param inData If input="raw", provide rgData here (Noob and dye-bias corrected signals produced by using `correct_noob_dye`). Otherwise, provide beta or M-value matrix here.
#' @param input Input data types. One of "raw", "beta", and "mval". If input is "beta" or "mval", please use probes as rows and samples as columns.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param vcfName VCF file name. Only effective when vcf=TRUE.
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
callGeno_typeII <- function(inData, input="raw", plotBeta=FALSE, vcf=FALSE, vcfName="genotypes.typeII_ccs_probe.vcf", 
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
  
  # filter probes based on peak density and positions.
  if(train){
    mod <- getMod(beta, cpu=cpu)
    beta <- beta[dplyr::filter(mod, h_0.1==TRUE, loc_pass==TRUE, nmod==3)$CpG,]
  }else{
    mod <- dplyr::filter(probeInfo_typeII, h_0.1==TRUE, loc_pass==TRUE, nmod==3) %>% 
      dplyr::select(SNP, CpG, h_0.1, loc_pass, nmod, loc0, loc1, loc2)
    rownames(mod) <- mod$CpG
    beta <- beta[mod$CpG,]
  }
  
  # calculate RAI
  mod$pM <- sapply(2 * mod$loc1, function(x) min(x, 1))
  RAI <- 1 - ( beta / matrix(rep(mod[rownames(beta), "pM", drop=TRUE], ncol(beta)), nrow=nrow(beta)) )
  RAI[RAI < 0.01] <- 0.01 # if set to 0 or 1, it will fail in fitting beta distribution as GP will be NA in call_genotypes.R:: GP <<- GP / tmp
  RAI[RAI > 0.99] <- 0.99
  
  # call genotypes
  genotypes <- call_genotypes(RAI, pop=pop, type="typeII_ccs_probe", maxiter=maxiter, 
                                       bayesian=bayesian, platform=platform, verbose=verbose)
  if(plotBeta){plot_beta_distribution(genotypes, type="typeII_ccs_probe")}
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, 
                             GP_cutoff=GP_cutoff, outlier_cutoff=outlier_cutoff, missing_cutoff=missing_cutoff,
                             R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, 
                             MAF_cutoff=MAF_cutoff, HWE_cutoff=HWE_cutoff, 
                             type="typeII_ccs_probe", pop=pop, plotAF=FALSE, platform=platform)
  list(dosage=dosage, mod=mod, genotypes=genotypes)
}

