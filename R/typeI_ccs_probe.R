
#' Call genotypes for Type I CCS probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param vcfName VCF file name. Only effective when vcf=TRUE.
#' @param GP_cutoff Genotypes with the highest genotype probability < GP_cutoff will be treated as missing. Only non-missing genotypes will be used to calculate MAF, R2, and HWE p value.
#' @param R2_cutoff_up,R2_cutoff_down R-square cutoffs to filter variants (Variants with R-square > R2_cutoff_up or < R2_cutoff_down should be removed). Note that for VCF output, variants with R-square outside this range will be marked in the `FILTER` column. For the returned dosage matrix, variants with R-square outside this range will be removed.
#' @param MAF_cutoff A MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with MAF below the cutoff will be removed.
#' @param HWE_cutoff HWE p value cutoff to filter variants. Note that for VCF output, variants with HWE p value below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with HWE p value below the cutoff will be removed.
#' @param missing_cutoff Missing rate cutoff to filter variants. Note that for VCF output, variants with missing rate above the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with missing rate above the cutoff will be removed.
#' @param train If TRUE, will fit the distribution of RAI (Ratio of Alternative allele Intensity) and filter probes by number of peaks. If FALSE, will use predefined probe list.
#' @param cpu Number of CPU. Only effective when train=TRUE.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL. Only probes with MAF of matching population > 0.01 will be kept. Only effective when train=TRUE.
#' @param bayesian Use the Bayesian approach to calculate posterior genotype probabilities.
#' @param platform EPIC or 450K.
#' @param verbose Verbose mode: 0/1/2.
#' @return A list containing
#' \item{dosage}{A matrix of genotype calls. Variants with R2 or MAF beyond the cutoffs are removed. Genotypes with the highest genotype probability < GP_cutoff will be marked as NA.}
#' \item{genotypes}{A list containing RAI, shapes of the mixed beta distributions, prior probabilities that the RAI values belong to one of the three genotypes, proportion of RAI values being outlier (U), and genotype probability (GP).}
#' @export
callGeno_typeI <- function(rgData, plotBeta=FALSE, vcf=FALSE, vcfName="genotypes.typeI_ccs_probe.vcf", 
                           GP_cutoff=0.9, R2_cutoff_up=1.1, R2_cutoff_down=0.75, MAF_cutoff=0.01, HWE_cutoff=1e-6, missing_cutoff=0.1, 
                           train=TRUE, cpu=1, pop="EAS", bayesian=TRUE, platform="EPIC", verbose=1){
  if(!train & platform!="EPIC"){
    print("Error: train=FALSE only works with platform=EPIC.")
    return(NA)
  }
  if(!(pop %in% c("EAS", "AMR", "AFR", "EUR", "SAS", "ALL"))){
    stop("pop must be one of EAS, AMR, AFR, EUR, SAS, and ALL.")
  }
  tag_af <- paste0(pop, "_AF")
  if(platform=="EPIC"){
    data(probeInfo_typeI)
  }else{
    data(probeInfo_typeI_450K); probeInfo_typeI <- probeInfo_typeI_450K
  }
  
  # calculate RAI
  df <- dplyr::filter(probeInfo_typeI, .data[["CpG"]] %in% rownames(rgData[["AR"]]) & .data[[tag_af]]>0.01 & .data[[tag_af]]<0.99)
  dR <- dplyr::filter(df, Color=="Red")
  dG <- dplyr::filter(df, Color=="Grn")
  dR_AR <- rgData[["AR"]][dR$CpG,] # Red channel, ib
  dR_BR <- rgData[["BR"]][dR$CpG,] # Red channel, ib
  dR_AG <- rgData[["AG"]][dR$CpG,] # Red channel, oob
  dR_BG <- rgData[["BG"]][dR$CpG,] # Red channel, oob
  dG_AR <- rgData[["AR"]][dG$CpG,] # Grn channel, oob
  dG_BR <- rgData[["BR"]][dG$CpG,] # Grn channel, oob
  dG_AG <- rgData[["AG"]][dG$CpG,] # Grn channel, ib
  dG_BG <- rgData[["BG"]][dG$CpG,] # Grn channel, ib
  RAI <- rbind(
    pmax(dG_AR + dG_BR, 1) / pmax(dG_AR + dG_BR + dG_AG + dG_BG, 2),
    pmax(dR_AG + dR_BG, 1) / pmax(dR_AG + dR_BG + dR_AR + dR_BR, 2)
  )
  
  # filter probes based on peak density and positions.
  if(train){
    mod <- getMod(RAI, cpu=cpu)
    RAI <- RAI[dplyr::filter(mod, h_0.1==TRUE, loc_pass==TRUE)$CpG,]
  }else{
    mod <- dplyr::filter(probeInfo_typeI, h_0.1==TRUE, loc_pass==TRUE) %>%
      dplyr::select(SNP, CpG, h_0.1, loc_pass, nmod, loc0, loc1, loc2)
    rownames(mod) <- mod$CpG
    RAI <- RAI[mod$CpG,]
  }
  
  # call genotypes
  genotypes <- call_genotypes_bayesian(RAI, pop=pop, type="typeI_ccs_probe", maxiter=50, 
                                       bayesian=bayesian, platform=platform, verbose=verbose)
  if(plotBeta){plot_beta_distribution(genotypes, type="typeI_ccs_probe")}
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, GP_cutoff=GP_cutoff,
                             R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, 
                             MAF_cutoff=MAF_cutoff, HWE_cutoff=HWE_cutoff, missing_cutoff=missing_cutoff, 
                             type="typeI_ccs_probe", pop=pop, plotAF=FALSE, platform=platform)
  list(dosage=dosage, mod=mod, genotypes=genotypes)
}

