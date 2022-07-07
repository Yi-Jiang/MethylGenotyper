
#' Call genotypes for Type-II probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param input Input data types. One of "raw", "beta", and "mval". If input is "beta" or "mval", please use probes as rows and samples as columns.
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
callGeno_typeII <- function(rgData, input="raw", plotBeta=FALSE, vcf=FALSE, vcfName="genotypes.typeII_probe.vcf", R2_cutoff_up=1.1, R2_cutoff_down=0.7, MAF_cutoff=0.01, train=FALSE, cpu=1, pop="ALL"){
  # calculate RAI
  if(input=="raw"){
    RAI <- getRAI_typeII(rgData, pop=pop)
  }else if(input=="beta"){
    data(probeInfo_typeII)
    tag_af <- paste0(pop, "_AF")
    beta <- rgData[rownames(rgData) %in% probeInfo_typeII[probeInfo_typeII[,tag_af]>0.01 & probeInfo_typeII[,tag_af]<0.99, "CpG"], ]
    RAI <- 1 - beta
  }else if(input=="mval"){
    data(probeInfo_typeII)
    tag_af <- paste0(pop, "_AF")
    beta <- mval2beta(rgData[rownames(rgData) %in% probeInfo_typeII[probeInfo_typeII[,tag_af]>0.01 & probeInfo_typeII[,tag_af]<0.99, "CpG"], ])
    RAI <- 1 - beta
  }else{
    print("Error: Input data type must be one of raw, beta, and mval.")
    return(NA)
  }
  
  # filter probes based on peak density and positions.
  if(train){
    mod <- getMod(RAI, cpu=cpu)
    RAI <- RAI[filter(mod, h_0.1==TRUE, loc_pass==TRUE)$Name,]
  }else{
    RAI <- RAI[rownames(RAI) %in% filter(probeInfo_typeII, h_0.1==TRUE, loc_pass==TRUE)$CpG,]
  }
  
  # call genotypes
  genotypes = callGeno(RAI, learn=TRUE)
  if(plotBeta){
    plot_beta_distribution(genotypes, type="typeII_probe")
  }
  dosage <- format_genotypes(genotypes, vcf=vcf, vcfName=vcfName, R2_cutoff_up=R2_cutoff_up, R2_cutoff_down=R2_cutoff_down, MAF_cutoff=MAF_cutoff, type="typeII_probe")
  dosage
}


#' Get RAI (Ratio of Alternative allele Intensity) for Type-II probes
#' 
#' @param rgData Noob and dye-bias corrected signals produced by using `correct_noob_dye`.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL. Only probes with MAF of matching population > 0.01 will be kept. Only effective when train=TRUE.
#' @return RAI (Ratio of Alternative allele Intensity).
#' @export
getRAI_typeII = function(rgData, pop="ALL"){
  if(!(pop %in% c("EAS", "AMR", "AFR", "EUR", "SAS", "ALL"))){
    stop("pop must be one of EAS, AMR, AFR, EUR, SAS, and ALL.")
  }
  tag_af <- paste0(pop, "_AF")
  data(probeInfo_typeII)
  cg <- rownames(rgData[["AR"]])
  cg_IIR <- cg[cg %in% probeInfo_typeII[probeInfo_typeII[,"Group"]=="IIR" & probeInfo_typeII[,tag_af]>0.01 & probeInfo_typeII[,tag_af]<0.99, "CpG"]] # Type II, Alt allele match Red
  #cg_IIG <- cg[cg %in% probeInfo_typeII[probeInfo_typeII[,"Group"]=="IIG" & probeInfo_typeII[,tag_af]>0.01 & probeInfo_typeII[,tag_af]<0.99, "CpG"]] # Type II, Alt allele match Grn
  RAI <- rbind(
    rgData[["AR"]][cg_IIR,] / (rgData[["AG"]][cg_IIR,] + rgData[["AR"]][cg_IIR,])
    #rgData[["AG"]][cg_IIG,] / (rgData[["AG"]][cg_IIG,] + rgData[["AR"]][cg_IIG,])
  )
  RAI
}

