
#' Call genotypes for SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param input Input data types. One of "raw", "beta", and "mval". If input is "beta" or "mval", please use probes as rows and samples as columns.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param vcfName VCF file name. Only effective when vcf=TRUE.
#' @param GP_cutoff When calculating missing rate, genotypes with the highest genotype probability < GP_cutoff will be treated as missing.
#' @param outlier_cutoff Genotypers with probability of being outliers > outlier_cutoff will be set as missing.
#' @param missing_cutoff Missing rate cutoff to filter variants. Note that for VCF output, variants with missing rate above the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with missing rate above the cutoff will be removed.
#' @param R2_cutoff_up,R2_cutoff_down R-square cutoffs to filter variants (Variants with R-square > R2_cutoff_up or < R2_cutoff_down should be removed). Note that for VCF output, variants with R-square outside this range will be marked in the `FILTER` column. For the returned dosage matrix, variants with R-square outside this range will be removed.
#' @param MAF_cutoff A MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with MAF below the cutoff will be removed.
#' @param HWE_cutoff HWE p value cutoff to filter variants. Note that for VCF output, variants with HWE p value below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with HWE p value below the cutoff will be removed.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param bayesian Use the Bayesian approach to calculate posterior genotype probabilities.
#' @param platform EPIC or 450K.
#' @param verbose Verbose mode: 0/1/2.
#' @return A list containing
#' \item{dosage}{A matrix of genotype calls. Variants with R2s, HWE p values, MAFs, or missing rates beyond the cutoffs are removed.}
#' \item{genotypes}{A list containing RAI, shapes of the mixed beta distributions, prior probabilities that the RAI values belong to one of the three genotypes, proportion of RAI values being outlier (U), and genotype probability (GP).}
#' @export
callGeno_snp <- function(rgData, input="raw", plotBeta=FALSE, vcf=FALSE, vcfName="genotypes.snp_probe.vcf", 
                         GP_cutoff=0.9, outlier_cutoff=0.2, missing_cutoff=0.1, 
                         R2_cutoff_up=1.1, R2_cutoff_down=0.75, MAF_cutoff=0.01, HWE_cutoff=1e-6, 
                         pop="EAS", bayesian=TRUE, platform="EPIC", verbose=1){
  if(input=="raw"){
    RAI <- getRAI_snp(rgData, platform=platform)
  }else if(input=="beta"){
    if(platform=="EPIC"){
      data(probeInfo_snp)
    }else{
      data(probeInfo_snp_450K); probeInfo_snp <- probeInfo_snp_450K
    }
    beta <- rgData[rownames(rgData) %in% probeInfo_snp$CpG, ]
    if(nrow(beta)==0){print("No SNP probes found. Exit!"); return(NA)}
    RAI_typeI <- 1 - beta[rownames(beta) %in% probeInfo_snp[probeInfo_snp$Group %in% c("IAG", "IAR", "IIR"), "CpG"], ] # alternative alleles match unmethylated probes.
    RAI_typeII <- beta[rownames(beta) %in% probeInfo_snp[probeInfo_snp$Group %in% c("IBG", "IBR", "IIG"), "CpG"], ]
    RAI <- rbind(RAI_typeI, RAI_typeII)
  }else if(input=="mval"){
    if(platform=="EPIC"){
      data(probeInfo_snp)
    }else{
      data(probeInfo_snp_450K); probeInfo_snp <- probeInfo_snp_450K
    }
    beta <- mval2beta(rgData[rownames(rgData) %in% probeInfo_snp$CpG, ])
    if(nrow(beta)==0){print("No SNP probes found. Exit!"); return(NA)}
    RAI_typeI <- 1 - beta[rownames(beta) %in% probeInfo_snp[probeInfo_snp$Group %in% c("IAG", "IAR", "IIR"), "CpG"], ]
    RAI_typeII <- beta[rownames(beta) %in% probeInfo_snp[probeInfo_snp$Group %in% c("IBG", "IBR", "IIG"), "CpG"], ]
    RAI <- rbind(RAI_typeI, RAI_typeII)
  }else{
    print("Error: Input data type must be one of raw, beta, and mval.")
    return(NA)
  }
  genotypes <- call_genotypes_bayesian(RAI, pop=pop, type="snp_probe", maxiter=50, bayesian=bayesian, platform=platform, verbose=verbose)
  if(plotBeta){plot_beta_distribution(genotypes, type="snp_probe")}
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, 
                             GP_cutoff=GP_cutoff, outlier_cutoff=outlier_cutoff, missing_cutoff=missing_cutoff, 
                             R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, 
                             MAF_cutoff=MAF_cutoff, HWE_cutoff=HWE_cutoff, 
                             type="snp_probe", pop=pop, plotAF=FALSE, platform=platform)
  list(dosage=dosage, genotypes=genotypes)
}

#' Convert M values to beta values
#' 
#' @param mval M value matrix.
#' @return Beta value matrix.
#' @export
mval2beta <- function(mval){
  mval_power <- 2 ^ mval
  beta <- mval_power / (1 + mval_power)
  beta
}

#' Get RAI (Ratio of Alternative allele Intensity) for SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param platform EPIC or 450K.
#' @return RAI (Ratio of Alternative allele Intensity).
#' @export
getRAI_snp <- function(rgData, platform="EPIC"){
  if(platform=="EPIC"){
    data(probeInfo_snp)
  }else{
    data(probeInfo_snp_450K); probeInfo_snp <- probeInfo_snp_450K
  }
  cg <- rownames(rgData[["AR"]])
  cg_IAR <- cg[cg %in% probeInfo_snp[probeInfo_snp$Group=="IAR", "CpG"]] # Type I, Red channel, Alt allele match probeA
  cg_IBR <- cg[cg %in% probeInfo_snp[probeInfo_snp$Group=="IBR", "CpG"]] # Type I, Red channel, Alt allele match probeB
  cg_IAG <- cg[cg %in% probeInfo_snp[probeInfo_snp$Group=="IAG", "CpG"]] # Type I, Grn channel, Alt allele match probeA
  cg_IBG <- cg[cg %in% probeInfo_snp[probeInfo_snp$Group=="IBG", "CpG"]] # Type I, Grn channel, Alt allele match probeB
  cg_IIR <- cg[cg %in% probeInfo_snp[probeInfo_snp$Group=="IIR", "CpG"]] # Type II, Alt allele match Red
  cg_IIG <- cg[cg %in% probeInfo_snp[probeInfo_snp$Group=="IIG", "CpG"]] # Type II, Alt allele match Grn
  RAI <- rbind(
    rgData[["AR"]][cg_IAR,] / (rgData[["AR"]][cg_IAR,] + rgData[["BR"]][cg_IAR,]),
    rgData[["BR"]][cg_IBR,] / (rgData[["AR"]][cg_IBR,] + rgData[["BR"]][cg_IBR,]),
    rgData[["AG"]][cg_IAG,] / (rgData[["AG"]][cg_IAG,] + rgData[["BG"]][cg_IAG,]),
    rgData[["BG"]][cg_IBG,] / (rgData[["AG"]][cg_IBG,] + rgData[["BG"]][cg_IBG,]),
    rgData[["AR"]][cg_IIR,] / (rgData[["AG"]][cg_IIR,] + rgData[["AR"]][cg_IIR,]),
    rgData[["AG"]][cg_IIG,] / (rgData[["AG"]][cg_IIG,] + rgData[["AR"]][cg_IIG,])
  )
  RAI
}

