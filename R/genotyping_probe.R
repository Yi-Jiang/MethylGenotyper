
#' Call genotypes for genotyping probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param R2_cutoff An R-square cutoff to filter variants. Note that for VCF output, variants with R-square below the cutoff will be marked in the `FILTER` column. For the returned matrix, variants with R-square below the cutoff will be removed.
#' @param MAF_cutoff An MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned matrix, variants with MAF below the cutoff will be removed.
#' @return A matrix of genotype calls.
#' @export
callGeno_genotyping <- function(rgData, plotBeta=FALSE, vcf=FALSE, R2_cutoff=0.7, MAF_cutoff=0.01){
  AB_geno <- getAB_genotyping(rgData)
  genotypes = ewastools::call_genotypes(AB_geno, learn=TRUE)
  if(plotBeta){
    plot_beta_distribution(genotypes, type="genotyping_probe")
  }
  dosage <- format_genotypes(genotypes, vcf=vcf, R2_cutoff=R2_cutoff, MAF_cutoff=MAF_cutoff, type="genotyping_probe")
  dosage
}


#' Get allelic balances for genotyping probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @return Allelic balances.
#' @export
getAB_genotyping = function(rgData){
  data(probeInfo_geno)
  cg <- rownames(rgData[["AR"]])
  cg_IAR <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IAR", "CpG"]] # Type I, Red channel, Alt allele match probeA
  cg_IBR <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IBR", "CpG"]] # Type I, Red channel, Alt allele match probeB
  cg_IAG <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IAG", "CpG"]] # Type I, Grn channel, Alt allele match probeA
  cg_IBG <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IBG", "CpG"]] # Type I, Grn channel, Alt allele match probeB
  cg_IIR <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IIR", "CpG"]] # Type II, Alt allele match Red
  cg_IIG <- cg[cg %in% probeInfo_geno[probeInfo_geno$Group=="IIG", "CpG"]] # Type II, Alt allele match Grn
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

