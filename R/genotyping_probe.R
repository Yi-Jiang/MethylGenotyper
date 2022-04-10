
#' Call genotypes for genotyping probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @return A matrix of genotype calls.
callGeno_genotyping <- function(rgData, plotBeta=FALSE, vcf=FALSE){
  AB_geno <- getAB_genotyping(rgData)
  genotypes = ewastools::call_genotypes(AB_geno, learn=TRUE)
  if(plotBeta){
    plot_beta_distribution(genotypes, type="genotyping_probe")
  }
  dosage <- format_genotypes(genotypes, vcf=vcf, type="genotyping_probe")
  dosage
}


#' Get allelic balances for genotyping probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @return Allelic balances.
getAB_genotyping = function(rgData){
  data(probeInfo_geno)
  cg <- rownames(rgData[["AR"]])
  cg_IAR <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IAR", "CpG"]] # Type I, Red channel, Alt allele match probeA
  cg_IBR <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IBR", "CpG"]] # Type I, Red channel, Alt allele match probeB
  cg_IAG <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IAG", "CpG"]] # Type I, Grn channel, Alt allele match probeA
  cg_IBG <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IBG", "CpG"]] # Type I, Grn channel, Alt allele match probeB
  cg_IIR <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IIR", "CpG"]] # Type II, Alt allele match Red
  cg_IIG <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IIG", "CpG"]] # Type II, Alt allele match Grn
  # cg_IAR <- cg[cg %in% grpGeno[["IAR"]]] # Type I, Red channel, Alt allele match probeA
  # cg_IBR <- cg[cg %in% grpGeno[["IBR"]]] # Type I, Red channel, Alt allele match probeB
  # cg_IAG <- cg[cg %in% grpGeno[["IAG"]]] # Type I, Grn channel, Alt allele match probeA
  # cg_IBG <- cg[cg %in% grpGeno[["IBG"]]] # Type I, Grn channel, Alt allele match probeB
  # cg_IIR <- cg[cg %in% grpGeno[["IIR"]]] # Type II, Alt allele match Red
  # cg_IIG <- cg[cg %in% grpGeno[["IIG"]]] # Type II, Alt allele match Grn
  AB_geno <- rbind(
    rgData[["AR"]][cg_IAR,] / (rgData[["AR"]][cg_IAR,] + rgData[["BR"]][cg_IAR,]),
    rgData[["BR"]][cg_IBR,] / (rgData[["AR"]][cg_IBR,] + rgData[["BR"]][cg_IBR,]),
    rgData[["AG"]][cg_IAG,] / (rgData[["AG"]][cg_IAG,] + rgData[["BG"]][cg_IAG,]),
    rgData[["BG"]][cg_IBG,] / (rgData[["AG"]][cg_IBG,] + rgData[["BG"]][cg_IBG,]),
    rgData[["AR"]][cg_IIR,] / (rgData[["AG"]][cg_IIR,] + rgData[["AR"]][cg_IIR,]),
    rgData[["AG"]][cg_IIG,] / (rgData[["AG"]][cg_IIG,] + rgData[["AR"]][cg_IIG,])
  )
  AB_geno
}

