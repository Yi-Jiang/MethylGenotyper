
#' Plot beta distributions for reference homozygous, heterozygous, and alternative homozygous
#' 
#' @param genotypes Genotype calls.
#' @param type One of "snp_probe", "typeI_ccs_probe", and "typeII_ccs_probe".
#' @export
plot_beta_distribution <- function(genotypes, type){
  shapes = as.matrix(genotypes$shapes[, c("shape1", "shape2")])
  pdf(paste0("beta_distribution.", type, ".pdf"), width=5, height=5)
  plot(
    x=seq(0, 1, 0.01), 
    y=dbeta(seq(0, 1, 0.01), shapes[1, "shape1"], shapes[1, "shape2"]), 
    type="l", main="Fitting model", xlab="Beta", ylab="Density", 
    cex.lab=1.8, cex.axis=1.5, xlim=c(0,1), ylim=c(0,30)
  )
  lines(x=seq(0, 1, 0.01), y=dbeta(seq(0, 1, 0.01), shapes[2, "shape1"], shapes[2, "shape2"]))
  lines(x=seq(0, 1, 0.01), y=dbeta(seq(0, 1, 0.01), shapes[3, "shape1"], shapes[3, "shape2"]))
  dev.off()
}

#' Constrain R2
#' 
#' R2 is calculated by var(G)/2p(1-p), where G is dosage genotype and p is allele frequency. Variants with 1<R2<=1.1 are constrained to 1. Variants with R2>1.1 (marked as .) are recommended to remove.
#' 
#' @param R2 R-square
#' @return Constrained R2
#' @export
constrain_R2 <- function(R2){
  if(R2>1){
    if(R2>1.1){
      R2 = "."
    }else{
      R2 = 1
    }
  }else{
    R2 <- sapply(R2, function(x) sprintf("%#.3g", x))
  }
  R2
}

#' Filter by R2
#' 
#' Variants with R2>1.1 (marked as .) or R2<0.75 are recommended to remove.
#' 
#' @param R2 R-square
#' @param R2_cutoff_up Variants with R-square greater than this cutoff should be removed.
#' @param R2_cutoff_down Variants with R-square less than this cutoff should be removed.
#' @return Whether the variant passed filtering.
#' @export
filter_by_R2 <- function(R2, R2_cutoff_up=1.1, R2_cutoff_down=0.75){
  if(R2>R2_cutoff_up){
    filter="R2_high"
  }else if(R2<R2_cutoff_down){
    filter="R2_low"
  }else{
    filter="PASS"
  }
  filter
}

#' Filter by AF
#' 
#' Variants with MAF<0.01 are recommended to remove.
#' 
#' @param AF Allele frequency
#' @param MAF_cutoff An MAF (Minor allele frequency) cutoff to filter variants.
#' @return Whether the variant passed filtering.
#' @export
filter_by_AF <- function(AF, MAF_cutoff=0.01){
  if(AF < MAF_cutoff | AF > (1-MAF_cutoff)){
    filter="MAF"
  }else{
    filter="PASS"
  }
  filter
}

#' Calculate Hardy-Weinberg Equilibrium (HWE) p value
#' 
#' @param hardgeno A matrix of hard genotypes, with each row indicates a SNP and each column indicates a sample.
#' @return HWE p values.
#' @export
getHWE <- function(hardgeno){
  if(is.null(rownames(hardgeno))){
    print("CAUTION: getHWE: It's recommended to provide rownames of the input matrix.")
  }
  hardgeno_sum <- t(apply(hardgeno, 1, function(x) table(x)[c("0/0", "0/1", "1/1")]))
  hardgeno_sum[is.na(hardgeno_sum)] <- 0
  colnames(hardgeno_sum) <- c("MM","MN","NN")
  hardgeno_sum_1 <- hardgeno_sum[rowSums(hardgeno_sum)!=0,]
  hwe <- suppressWarnings(HWChisqMat(hardgeno_sum_1, cc=0, verbose=FALSE))
  if(is.null(rownames(hardgeno))){
    hwe_p <- hwe$pvalvec
  }else{
    hwe_p <- rep(0, nrow(hardgeno_sum))
    names(hwe_p) <- rownames(hardgeno_sum)
    hwe_p[rownames(hardgeno_sum_1)] <- hwe$pvalvec
  }
  hwe_p
}

#' Filter by Hardy–Weinberg Equilibrium (HWE) p value
#' 
#' Variants with Hardy–Weinberg Equilibrium (HWE) p value < HWE_cutoff are recommended to remove.
#' 
#' @param hwe_p Hardy–Weinberg Equilibrium (HWE) p values
#' @param HWE_cutoff A HWE p value cutoff to filter variants.
#' @return Whether the variant passed filtering.
#' @export
filter_by_HWE <- function(hwe_p, HWE_cutoff=1e-6){
  if(hwe_p < HWE_cutoff){
    filter="HWE"
  }else{
    filter="PASS"
  }
  filter
}

#' Filter by missing rate
#' 
#' Variants with missing rate > missing_cutoff are recommended to remove.
#' 
#' @param F_MISSING Fraction of missing genotypes
#' @param missing_cutoff Missing rate cutoff to filter variants.
#' @return Whether the variant passed filtering.
#' @export
filter_by_missing <- function(F_MISSING, missing_cutoff=0.1){
  if(F_MISSING > missing_cutoff){
    filter="Missing"
  }else{
    filter="PASS"
  }
  filter
}

#' Get hard genotypes from genotype probabilities
#' 
#' Hard genotype is defined as the genotype with highest genotype probability.
#' 
#' @param AA Genotype probability of AA or 0/0.
#' @param AB Genotype probability of AB or 0/1.
#' @param BB Genotype probability of BB or 1/1.
#' @return A matrix of hard genotypes.
#' @export
dosage2hard <- function(AA, AB, BB){
  pmax0 <- pmax(AA, AB, BB)
  hardgeno <- matrix(NA, nrow=nrow(AA), ncol=ncol(AA), dimnames=dimnames(AA))
  hardgeno[pmax0==AA] <- "0/0"
  hardgeno[pmax0==AB] <- "0/1"
  hardgeno[pmax0==BB] <- "1/1"
  #hardgeno[pmax0 < GP_cutoff] <- NA_real_
  hardgeno
}

#' Format genotype calls
#' 
#' @param genotypes Genotype calls.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param vcfName VCF file name. Only effective when vcf=TRUE.
#' @param GP_cutoff When calculating missing rate, genotypes with the highest genotype probability < GP_cutoff will be treated as missing.
#' @param outlier_cutoff "max" or a number ranging from 0 to 1. If outlier_cutoff="max", genotypes with outlier probability larger than all of the three genotype probabilities will be set as missing. If outlier_cutoff is a number, genotypes with outlier probability > outlier_cutoff will be set as missing.
#' @param missing_cutoff Missing rate cutoff to filter variants. Note that for VCF output, variants with missing rate above the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with missing rate above the cutoff will be removed.
#' @param R2_cutoff_up,R2_cutoff_down R-square cutoffs to filter variants (Variants with R-square > R2_cutoff_up or < R2_cutoff_down should be removed). Note that for VCF output, variants with R-square outside this range will be marked in the `FILTER` column. For the returned dosage matrix, variants with R-square outside this range will be removed.
#' @param MAF_cutoff MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with MAF below the cutoff will be removed.
#' @param HWE_cutoff HWE p value cutoff to filter variants. Note that for VCF output, variants with HWE p value below the cutoff will be marked in the `FILTER` column. For the returned dosage matrix, variants with HWE p value below the cutoff will be removed.
#' @param pop Population to be used to extract AFs. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @param plotAF To plot the distribution of AFs in 1KGP and input data.
#' @param platform EPIC or 450K.
#' @return A matrix of genotype calls. Variants with R2s, HWE p values, MAFs, or missing rates beyond the cutoffs are removed.
#' @export
format_genotypes <- function(genotypes, vcf=FALSE, vcfName, GP_cutoff=0.9, outlier_cutoff="max", missing_cutoff=0.1, 
                             R2_cutoff_up=1.1, R2_cutoff_down=0.75, MAF_cutoff=0.01, HWE_cutoff=1e-6, 
                             pop="ALL", type, plotAF=FALSE, platform="EPIC"){
  print(paste(Sys.time(), "Calculating AF, R2, and HWE."))
  dosage <- genotypes$GP$pAB + 2 * genotypes$GP$pBB
  probes <- rownames(dosage)
  maxGP <- pmax(genotypes$GP$pAA, genotypes$GP$pAB, genotypes$GP$pBB, na.rm=TRUE)
  #dosage[maxGP < GP_cutoff] <- NA_real_
  if(outlier_cutoff=="max"){
    dosage[genotypes$outliers > maxGP * (1 - genotypes$outliers)] <- NA_real_
  }else{
    dosage[genotypes$outliers > outlier_cutoff] <- NA_real_
  }
  missing <- rowSums(maxGP < GP_cutoff | is.na(dosage)) / ncol(maxGP)
  AF <- rowMeans(dosage, na.rm=T) / 2
  AF[is.na(AF)] <- 0
  R2 <- apply(dosage, 1, function(x) var(x, na.rm=T)) / (2 * AF * (1 - AF))
  R2[is.na(R2)] <- 0
  R2_constrained <- sapply(R2, constrain_R2)
  hardgeno <- dosage2hard(genotypes$GP$pAA, genotypes$GP$pAB, genotypes$GP$pBB)
  hardgeno[is.na(dosage)] <- NA_real_
  hwe_p <- getHWE(hardgeno)
  filter_AF <- sapply(AF, function(x) filter_by_AF(x, MAF_cutoff))
  filter_R2 <- sapply(R2, function(x) filter_by_R2(x, R2_cutoff_up, R2_cutoff_down))
  filter_HWE <- sapply(hwe_p, function(x) filter_by_HWE(x, HWE_cutoff))
  filter_missing <- sapply(missing, function(x) filter_by_missing(x, missing_cutoff))
  filter <- paste(filter_AF, filter_R2, filter_HWE, filter_missing, sep=";")
  filter[filter=="PASS;PASS;PASS;PASS"] <- "PASS"
  filter[filter!="PASS;PASS;PASS;PASS"] <- gsub(";PASS", "", gsub("PASS;", "", filter[filter!="PASS;PASS;PASS;PASS"]))
  AF <- sapply(AF, function(x) sprintf("%#.3g", x))

  ## Write into a VCF file
  if(vcf){
    ## format
    hardgeno[is.na(hardgeno)] <- "./."
    dosage_chr <- round(dosage, 2)
    dosage_chr[is.na(dosage_chr)] <- "."
    genotypes$GP$pAA <- round(genotypes$GP$pAA, 2)
    genotypes$GP$pAB <- round(genotypes$GP$pAB, 2)
    genotypes$GP$pBB <- round(genotypes$GP$pBB, 2)

    ## Genotype
    print(paste(Sys.time(), "Preparing and writing VCF file."))
    geno <- matrix(nrow=nrow(dosage), ncol=ncol(dosage))
    colnames(geno) <- colnames(dosage)
    for(i in 1:nrow(geno)){
      for(j in 1:ncol(geno)){
        geno[i,j] <- paste0(
          hardgeno[i,j], ":",
          dosage_chr[i,j], ":", 
          round(genotypes$RAI[i,j], 2), ":", 
          genotypes$GP$pAA[i,j], ",", 
          genotypes$GP$pAB[i,j], ",", 
          genotypes$GP$pBB[i,j])
      }
    }
    geno <- as.data.frame(cbind(CpG = rownames(dosage), geno))
    if(type=="snp_probe"){
      if(platform=="EPIC"){
        data(probeInfo_snp); probeInfo <- probeInfo_snp
      }else{
        data(probeInfo_snp_450K); probeInfo <- probeInfo_snp_450K
      }
    }else if(type=="typeI_ccs_probe"){
      if(platform=="EPIC"){
        data(probeInfo_typeI); probeInfo <- probeInfo_typeI
      }else{
        data(probeInfo_typeI_450K); probeInfo <- probeInfo_typeI_450K
      }
    }else if(type=="typeII_ccs_probe"){
      if(platform=="EPIC"){
        data(probeInfo_typeII); probeInfo <- probeInfo_typeII
      }else{
        data(probeInfo_typeII_450K); probeInfo <- probeInfo_typeII_450K
      }
    }else{
      print("Error: misspecified probe type!")
      return
    }
    rownames(probeInfo) <- probeInfo$CpG
    vcf <- cbind(
      probeInfo[probes, 1:6], 
      tibble(QUAL=".", FILTER=filter, 
             INFO=paste0("AF=", AF, ";R2=", R2_constrained, ";HWE=", sprintf("%.3g", hwe_p), ";F_MISSING=", sprintf("%.2g", missing)), 
             FORMAT="GT:DS:RAI:GP")
    ) %>% left_join(geno, by=c("CpG"))
    vcf <- vcf[, -6]
    vcf$Chr <- factor(vcf$Chr, levels=c(paste0("chr", 1:22)))
    vcf <- vcf[order(vcf$Chr, vcf$Pos),]
    
    ## Header
    samplelist <- paste(colnames(dosage), collapse="\t")
    header <- paste(
      "##fileformat=VCFv4.2",
      paste(paste0("##contig=<ID=chr", 1:22, ">"), collapse="\n"),
      paste0("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">"),
      paste0("##INFO=<ID=R2,Number=1,Type=Float,Description=\"R-square, encoded as var(G)/2p(1-p), where G is dosage genotype and p is allele frequency. Variants with 1<R2<=1.1 are constrained to 1. Variants with R2>1.1 (marked as .) are recommended to remove\">"),
      paste0("##INFO=<ID=HWE,Number=1,Type=Float,Description=\"Hardy-Weinberg Equilibrium p-value\">"),
      paste0("##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes, which is defined as genotypes with the highest genotype probability < ", GP_cutoff, "\">"),
      paste0("##FILTER=<ID=PASS,Description=\"All filters passed\">"),
      paste0("##FILTER=<ID=MAF,Description=\"MAF is below ", MAF_cutoff, "\">"),
      paste0("##FILTER=<ID=R2_low,Description=\"R2 is below ", R2_cutoff_down, "\">"),
      paste0("##FILTER=<ID=R2_high,Description=\"R2 is above ", R2_cutoff_up, "\">"),
      paste0("##FILTER=<ID=HWE,Description=\"Deviation from Hardy-Weinberg Equilibrium (HWE, p < 1E-6)\">"),
      paste0("##FILTER=<ID=Missing,Description=\"Fraction of missing genotypes is above ", missing_cutoff, "\">"),
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage\">",
      "##FORMAT=<ID=RAI,Number=1,Type=Float,Description=\"RAI (Ratio of Alternative allele Intensity)\">",
      "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probability of reference homozygous, heterozygous, and alternative homozygous\">",
      paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", samplelist),
      sep = "\n"
    )
    write.table(header, file=vcfName, sep="\t", row.names=F, quote=F, col.names=F)
    write.table(vcf, file=vcfName, sep="\t", row.names=F, quote=F, col.names=F, append=T)
  }
  
  ## Filter dosage
  dosage <- dosage[filter=="PASS",,drop=F]
  
  ## Plots
  probe2af <- get_AF(pop=pop, type=type, platform=platform)
  if(plotAF){plotAF_func(AF_input=AF[rownames(dosage)], AF_1KGP=probe2af[rownames(dosage)], pop=pop, type=type)}
  
  dosage
}

#' Plot the distribution of AFs in 1KGP and input data.
#'
#' @param AF_input A vector.
#' @param AF_1KGP A vector.
#' @param pop Population. One of EAS, AMR, AFR, EUR, SAS, and ALL.
#' @param type One of snp_probe, typeI_ccs_probe, and typeII_ccs_probe.
#' @export
plotAF_func <- function(AF_input, AF_1KGP, pop, type){
  probes <- intersect(names(AF_input), names(AF_1KGP))
  d <- data.frame(AF_input=AF_input[probes], AF_1KGP=AF_1KGP[probes])
  p <- ggplot(d) + 
    geom_point(aes(x=AF_input, y=AF_1KGP), alpha=.5) + 
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    theme_bw() + labs(x="AF (input)", y=paste0("AF (1KGP ", pop, ")"), title=type)
  ggsave(filename=paste0("AF.", type, ".pdf"), plot=p, width=5, height=5, units="in", scale=2)
}
