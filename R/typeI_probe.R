
#' Call genotypes for Type I CCS probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotRAI If TRUE, plot distribution of RAIs.
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
#' @param bw band width.
#' @param minDens A parameter for mode test. Minimum density for a valid peak.
#' @param maxProp_antimode A parameter for mode test. If the proportion of antimode density on mode density exceed this threshold, the two peaks will be merged.
#' @param bayesian Use the Bayesian approach to calculate posterior genotype probabilities.
#' @param platform EPIC or 450K.
#' @param verbose Verbose mode: 0/1/2.
#' @return A list containing
#' \item{dosage}{A matrix of genotype calls. Variants with R2s, HWE p values, MAFs, or missing rates beyond the cutoffs are removed.}
#' \item{genotypes}{A list containing RAI, shapes of the mixed beta distributions, prior probabilities that the RAI values belong to one of the three genotypes, proportion of RAI values being outlier (U), and genotype probability (GP).}
#' @export
callGeno_typeI <- function(rgData, plotRAI=FALSE, vcf=FALSE, vcfName="genotypes.typeI_ccs_probe.vcf", 
                           bw=0.04, minDens=0.01, maxProp_antimode=0.5,
                           GP_cutoff=0.9, outlier_cutoff="max", missing_cutoff=0.1, 
                           R2_cutoff_up=1.1, R2_cutoff_down=0.75, MAF_cutoff=0.01, HWE_cutoff=1e-6, 
                           train=TRUE, cpu=1, pop="EAS", bayesian=FALSE, platform="EPIC", verbose=1){
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
  
  # remove probes if they have common SNPs (MAF>0.01 in corresponding population) within 5bps
  probeInfo_typeI <- probeInfo_typeI[!is.na(probeInfo_typeI[,tag_af]),]
  
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
    mod <- getMod(RAI, bw=bw, minDens=minDens, maxProp_antimode=maxProp_antimode, cpu=cpu)
    RAI <- RAI[dplyr::filter(mod, loc_pass==TRUE, nmod>=3)$CpG,]
  }else{
    mod <- dplyr::filter(probeInfo_typeI, loc_pass==TRUE, nmod>=2) %>%
      dplyr::select(SNP, CpG, h_0.1, loc_pass, nmod, loc0, loc1, loc2)
    rownames(mod) <- mod$CpG
    RAI <- RAI[mod$CpG,]
  }
  
  # call genotypes
  genotypes <- call_genotypes(RAI, pop=pop, type="typeI_ccs_probe", maxiter=50, 
                                       bayesian=bayesian, platform=platform, verbose=verbose)
  if(plotRAI){plot_RAI_distribution(genotypes, type="typeI_ccs_probe")}
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, 
                             GP_cutoff=GP_cutoff, outlier_cutoff=outlier_cutoff, missing_cutoff=missing_cutoff,
                             R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, 
                             MAF_cutoff=MAF_cutoff, HWE_cutoff=HWE_cutoff, 
                             type="typeI_ccs_probe", pop=pop, plotAF=FALSE, platform=platform)
  list(dosage=dosage, mod=mod, genotypes=genotypes)
}

