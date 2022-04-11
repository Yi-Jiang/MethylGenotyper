
#' Call genotypes for Type-II SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param plotBeta If TRUE, plot beta distributions for reference homozygous, heterozygous, and alternative homozygous.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param R2_cutoff An R-square cutoff to filter variants. Note that for VCF output, variants with R-square below the cutoff will be marked in the `FILTER` column. For the returned matrix, variants with R-square below the cutoff will be removed.
#' @param MAF_cutoff An MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned matrix, variants with MAF below the cutoff will be removed.
#' @return A matrix of genotype calls.
#' @export
callGeno_typeII <- function(rgData, plotBeta=FALSE, vcf=FALSE, R2_cutoff=0.7, MAF_cutoff=0.01){
  AB_geno <- getAB_typeII(rgData)
  genotypes = ewastools::call_genotypes(AB_geno, learn=TRUE)
  if(plotBeta){
    plot_beta_distribution(genotypes, type="typeII_snp_probe")
  }
  dosage <- format_genotypes(genotypes, vcf=vcf, R2_cutoff=R2_cutoff, MAF_cutoff=MAF_cutoff, type="typeII_snp_probe")
  dosage
}


#' Get allelic balances for Type-II SNP probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @return Allelic balances.
#' @export
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

