
#' Plot beta distributions for reference homozygous, heterozygous, and alternative homozygous
#' 
#' @param genotypes Genotype calls produced by using `ewastools::call_genotypes`.
#' @param type One of "genotyping_probe", "ccs_snp_probe", and "typeII_snp_probe".
plot_beta_distribution <- function(genotypes, type){
  pdf(paste0("beta_distribution.", type, ".pdf"), width=500, height=500)
  plot(
    x=seq(0, 1, 0.01), 
    y=dbeta(seq(0, 1, 0.01), genotypes$par$shapes1[1], genotypes$par$shapes2[1]), 
    type="l", main="Fitting model", xlab="Beta", ylab="Density", 
    cex.lab=1.8, cex.axis=1.5, xlim=c(0,1), ylim=c(0,30)
  )
  lines(x=seq(0, 1, 0.01), y=dbeta(seq(0, 1, 0.01), genotypes$par$shapes1[2], genotypes$par$shapes2[2]))
  lines(x=seq(0, 1, 0.01), y=dbeta(seq(0, 1, 0.01), genotypes$par$shapes1[3], genotypes$par$shapes2[3]))
  dev.off()
}

#' Constrain R2
#' 
#' R2 is calculated by var(G)/2p(1-p), where G is dosage genotype and p is allele frequency. Variants with 1<R2<=1.1 are constrained to 1. Variants with R2>1.1 (marked as .) are recommended to remove.
#' 
#' @param R2 R-square
#' @return Constrained R2
constrain_R2 <- function(R2){
  if(R2>1){
    if(R2>1.1){
      R2 = "."
    }else{
      R2 = 1
    }
  }else{
    R2 = round(R2, 3)
  }
  R2
}

#' Filter by R2
#' 
#' Variants with R2>1.1 (marked as .) or R2<0.7 are recommended to remove.
#' 
#' @param R2 R-square
#' @return Whether the variant passed filtering.
filter_by_R2 <- function(R2){
  if(R2>1.1){
    filter="R2>1.1"
  }else if(R2<0.7){
    filter="R2<0.7"
  }else{
    filter="PASS"
  }
  filter
}

#' Format genotype calls produced by ewastools
#' 
#' @param genotypes Genotype calls produced by using `ewastools::call_genotypes`.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param type One of genotyping_probe, ccs_snp_probe, and typeII_snp_probe.
#' @return A matrix of genotype calls.
format_genotypes <- function(genotypes, vcf=FALSE, type){
  for(i in 1:3){
    colnames(genotypes$gamma[[i]]) <- colnames(genotypes$snps)
    rownames(genotypes$gamma[[i]]) <- rownames(genotypes$snps)
  }
  dosage <- genotypes$gamma[[2]] + 2 * genotypes$gamma[[3]]
  AF = rowMeans(dosage) / 2
  R2 <- apply(dosage, 1, var) / (2 * AF * (1 - AF))
  R2[is.na(R2)] <- 0
  R2_constrained <- sapply(R2, constrain_R2)
  filter <- sapply(R2, filter_by_R2)
  
  ## Write into a VCF file
  if(vcf){
    samplelist <- paste(colnames(dosage), collapse="\t")
    header <- paste(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
      "##INFO=<ID=R2,Number=1,Type=Float,Description=\"R-square, encoded as var(G)/2p(1-p), where G is dosage genotype and p is allele frequency. Variants with 1<R2<=1.1 are constrained to 1. Variants with R2>1.1 (marked as .) are recommended to remove.\">",
      "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage.\">",
      "##FORMAT=<ID=AB,Number=1,Type=Float,Description=\"Allelic balance.\">",
      "##FORMAT=<ID=GP,Number=1,Type=Float,Description=\"Genotype probability of reference homozygous, heterozygous, and alternative homozygous produced by using ewastools::call_genotypes.\">",
      paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", samplelist),
      sep = "\n"
    )
    geno <- matrix(nrow=nrow(dosage), ncol=ncol(dosage))
    colnames(geno) <- colnames(dosage)
    for(i in 1:nrow(geno)){
      for(j in 1:ncol(geno)){
        geno[i,j] <- paste0(
          round(dosage[i,j], 2), ":", 
          round(genotypes$snps[i,j], 2), ":", 
          round(genotypes$gamma[[1]][i,j], 2), ",", 
          round(genotypes$gamma[[2]][i,j], 2), ",", 
          round(genotypes$gamma[[3]][i,j], 2))
      }
    }
    geno <- as.data.frame(cbind(CpG = rownames(dosage), geno))
    if(type=="genotyping_probe"){
      data(probeInfo_geno); probeInfo = probeInfo_geno
    }else if(type=="ccs_snp_probe"){
      data(probeInfo_ccs); probeInfo = probeInfo_ccs
    }else if(type=="typeII_snp_probe"){
      data(probeInfo_typeII); probeInfo = probeInfo_typeII
    }else{
      print("Error: misspecified probe types!")
      return
    }
    vcf <- cbind(
      probeInfo[,1:6], 
      tibble(QUAL=".", FILTER=filter, INFO=paste0("AF=", round(AF, 3), ";R2=", R2_constrained), FORMAT="DS:AB:GP")
    ) %>% left_join(geno, by=c("CpG"))
    vcf <- vcf[, -6]
    write.table(header, file=paste0("genotypes.", type, ".vcf"), sep="\t", row.names=F, quote=F, col.names=F)
    write.table(vcf, file=paste0("genotypes.", type, ".vcf"), sep="\t", row.names=F, quote=F, col.names=F, append=T)
    dosage
  }
}
