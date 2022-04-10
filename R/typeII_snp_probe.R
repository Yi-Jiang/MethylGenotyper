
#' Call genotypes for Type-II SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @return A matrix of genotype calls.
callGeno_typeII <- function(rgData, plotBeta=FALSE, vcf=FALSE){
  AB_geno <- getAB_typeII(rgData)
  genotypes = ewastools::call_genotypes(AB_geno, learn=TRUE)
  if(plotBeta){
    plot_beta_distribution(genotypes, type="typeII_snp_probe")
  }
  dosage <- format_genotypes(genotypes, vcf=vcf, type="typeII_snp_probe")
  dosage
}


#' Get allelic balances for Type-II SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @return Allelic balances.
getAB_typeII = function(rgData){
  data(probeInfo_typeII)
  cg <- rownames(rgData[["AR"]])
  cg_IIR <- cg[cg %in% probeInfo_typeII[probeInfo_typeII$Group=="IIR", "CpG"]] # Type II, Alt allele match Red
  cg_IIG <- cg[cg %in% probeInfo_typeII[probeInfo_typeII$Group=="IIG", "CpG"]] # Type II, Alt allele match Grn
  AB_geno <- rbind(
    rgData[["AR"]][cg_IIR,] / (rgData[["AG"]][cg_IIR,] + rgData[["AR"]][cg_IIR,]),
    rgData[["AG"]][cg_IIG,] / (rgData[["AG"]][cg_IIG,] + rgData[["AR"]][cg_IIG,])
  )
  AB_geno
}

