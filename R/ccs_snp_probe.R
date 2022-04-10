
#' Call genotypes for CCS SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @return A matrix of genotype calls.
#' @export
callGeno_ccs <- function(rgData, plotBeta=FALSE, vcf=FALSE){
  AB_geno <- getAB_ccs(rgData)
  genotypes = ewastools::call_genotypes(AB_geno, learn=TRUE)
  if(plotBeta){
    plot_beta_distribution(genotypes, type="ccs_snp_probe")
  }
  dosage <- format_genotypes(genotypes, vcf=vcf, type="ccs_snp_probe")
  dosage
}


#' Get allelic balances for CCS SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @return Allelic balances.
#' @export
getAB_ccs = function(rgData){
  data(probeInfo_ccs)
  dR <- filter(probeInfo_ccs, Color=="Red")
  dG <- filter(probeInfo_ccs, Color=="Grn")
  dR_AR <- rgData[["AR"]][dR$CpG,] # Red channel, ib
  dR_BR <- rgData[["BR"]][dR$CpG,] # Red channel, ib
  dR_AG <- rgData[["AG"]][dR$CpG,] # Red channel, oob
  dR_BG <- rgData[["BG"]][dR$CpG,] # Red channel, oob
  dG_AR <- rgData[["AR"]][dG$CpG,] # Grn channel, oob
  dG_BR <- rgData[["BR"]][dG$CpG,] # Grn channel, oob
  dG_AG <- rgData[["AG"]][dG$CpG,] # Grn channel, ib
  dG_BG <- rgData[["BG"]][dG$CpG,] # Grn channel, ib
  AB_geno <- rbind(
    pmax(dG_AR + dG_BR, 1) / pmax(dG_AR + dG_BR + dG_AG + dG_BG, 2),
    pmax(dR_AG + dR_BG, 1) / pmax(dR_AG + dR_BG + dR_AR + dR_BR, 2)
  )
  AB_geno
}

