
#' Call genotypes for CCS probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param vcfName VCF file name. Only effective when vcf=TRUE.
#' @param R2_cutoff_up,R2_cutoff_down R-square cutoffs to filter variants (Variants with R-square > R2_cutoff_up or < R2_cutoff_down should be removed). Note that for VCF output, variants with R-square outside this range will be marked in the `FILTER` column. For the returned matrix, variants with R-square outside this range will be removed.
#' @param MAF_cutoff A MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned matrix, variants with MAF below the cutoff will be removed.
#' @param train If TRUE, will fit the distribution of RAI (Ratio of Alternative allele Intensity) and filter probes by number of peaks. If FALSE, will use predefined probe list.
#' @param cpu Number of CPU. Only effective when train=TRUE.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL. Only probes with MAF of matching population > 0.01 will be kept. Only effective when train=TRUE.
#' @return A matrix of genotype calls.
#' @export
callGeno_ccs <- function(rgData, plotBeta=FALSE, vcf=FALSE, vcfName="genotypes.ccs_snp_probe.vcf", R2_cutoff_up=1.1, R2_cutoff_down=0.7, MAF_cutoff=0.01, train=FALSE, cpu=1, pop="ALL"){
  RAI <- getRAI_ccs(rgData, pop=pop)
  if(train){
    mod <- getMod(RAI, cpu=cpu)
    RAI <- RAI[filter(mod, h_0.1==TRUE, loc_pass==TRUE)$Name,]
  }else{
    RAI <- RAI[rownames(RAI) %in% filter(probeInfo_ccs, h_0.1==TRUE, loc_pass==TRUE)$CpG,]
  }
  genotypes = callGeno(RAI, learn=TRUE)
  if(plotBeta){
    plot_beta_distribution(genotypes, type="ccs_snp_probe")
  }
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, MAF_cutoff=MAF_cutoff, type="ccs_snp_probe")
  dosage
}

#' Get RAI (Ratio of Alternative allele Intensity) for CCS probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL. Only probes with MAF of matching population > 0.01 will be kept. Only effective when train=TRUE.
#' @return RAI (Ratio of Alternative allele Intensity).
#' @export
getRAI_ccs = function(rgData, pop="ALL"){
  if(!(pop %in% c("EAS", "AMR", "AFR", "EUR", "SAS", "ALL"))){
    stop("pop must be one of EAS, AMR, AFR, EUR, SAS, and ALL.")
  }
  tag_af <- paste0(pop, "_AF")
  data(probeInfo_ccs)
  dR <- filter(probeInfo_ccs, Color=="Red", .data[[tag_af]]>0.01 & .data[[tag_af]]<0.99)
  dG <- filter(probeInfo_ccs, Color=="Grn", .data[[tag_af]]>0.01 & .data[[tag_af]]<0.99)
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
  RAI
}

