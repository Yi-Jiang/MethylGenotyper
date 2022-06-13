
#' Plot beta distributions for reference homozygous, heterozygous, and alternative homozygous
#' 
#' @param genotypes Genotype calls.
#' @param type One of "snp_probe", "ccs_snp_probe", and "typeII_probe".
#' @export
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
#' @export
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
#' @param R2_cutoff_up Variants with R-square greater than this cutoff should be removed.
#' @param R2_cutoff_down Variants with R-square less than this cutoff should be removed.
#' @return Whether the variant passed filtering.
#' @export
filter_by_R2 <- function(R2, R2_cutoff_up=1.1, R2_cutoff_down=0.7){
  if(R2>R2_cutoff_up){
    filter=paste0("R2>", R2_cutoff_up)
  }else if(R2<R2_cutoff_down){
    filter=paste0("R2<", R2_cutoff_down)
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
    filter=paste0("MAF<", MAF_cutoff)
  }else{
    filter="PASS"
  }
  filter
}

#' Calculate Hardy–Weinberg Equilibrium (HWE) p value
#' 
#' @param hardgeno A matrix of hard genotypes, with each row indicates a SNP and each column indicates a sample.
#' @return HWE p values.
#' @export
getHWE <- function(hardgeno){
  hardgeno_sum <- t(apply(hardgeno, 1, function(x) table(x)[c("0/0", "0/1", "1/1")]))
  hardgeno_sum[is.na(hardgeno_sum)] <- 0
  colnames(hardgeno_sum) <- c("MM","MN","NN")
  hardgeno_sum_1 <- hardgeno_sum[rowSums(hardgeno_sum)!=0,]
  hardgeno_sum_0 <- hardgeno_sum[rowSums(hardgeno_sum)==0,]
  hwe <- suppressWarnings(HWChisqMat(hardgeno_sum_1, cc=0, verbose=FALSE))
  hwe_p <- rep(0, nrow(hardgeno_sum))
  names(hwe_p) <- rownames(hardgeno_sum)
  hwe_p[rownames(hardgeno_sum_1)] <- hwe$pvalvec
  hwe_p
}

#' Filter by Hardy–Weinberg Equilibrium (HWE) p value
#' 
#' Variants with Hardy–Weinberg Equilibrium (HWE) p value < 1e-6 are recommended to remove.
#' 
#' @param hwe_p Hardy–Weinberg Equilibrium (HWE) p values
#' @param hwe_p_cutoff A HWE p value cutoff to filter variants.
#' @return Whether the variant passed filtering.
#' @export
filter_by_HWE <- function(hwe_p, hwe_p_cutoff){
  if(hwe_p < hwe_p_cutoff){
    filter=paste0("P(HWE)<", hwe_p_cutoff)
  }else{
    filter="PASS"
  }
  filter
}

#' Get hard genotypes from genotype probabilities
#' 
#' Hard genotype is defined as the genotype with highest genotype probability.
#' 
#' @param genotypes Genotype calls.
#' @return A matrix of hard genotypes.
#' @export
dosage2hard <- function(genotypes){
  hardgeno <- t(sapply(1:nrow(genotypes$snps), function(i){
    sapply(1:ncol(genotypes$snps), function(j){
      if(is.na(genotypes$gamma[[1]][i,j])){
        NA
      }else{
        names(which.max(c(
          "0/0"=genotypes$gamma[[1]][i,j], 
          "0/1"=genotypes$gamma[[2]][i,j], 
          "1/1"=genotypes$gamma[[3]][i,j]
        )))
      }
    })
  }))
  dimnames(hardgeno) <- dimnames(genotypes$snps)
  hardgeno
}

#' Format genotype calls produced by ewastools
#' 
#' @param genotypes Genotype calls.
#' @param vcf If TRUE, will write a VCF file in the current directory.
#' @param R2_cutoff_up,R2_cutoff_down R-square cutoffs to filter variants (Variants with R-square > R2_cutoff_up or < R2_cutoff_down should be removed). Note that for VCF output, variants with R-square outside this range will be marked in the `FILTER` column. For the returned matrix, variants with R-square outside this range will be removed.
#' @param MAF_cutoff An MAF cutoff to filter variants. Note that for VCF output, variants with MAF below the cutoff will be marked in the `FILTER` column. For the returned matrix, variants with MAF below the cutoff will be removed.
#' @param type One of snp_probe, ccs_snp_probe, and typeII_probe.
#' @return A matrix of genotype calls.
#' @export
format_genotypes <- function(genotypes, vcf=FALSE, R2_cutoff_up=1.1, R2_cutoff_down=0.7, MAF_cutoff=0.01, type){
  for(i in 1:3){
    colnames(genotypes$gamma[[i]]) <- colnames(genotypes$snps)
    rownames(genotypes$gamma[[i]]) <- rownames(genotypes$snps)
  }
  dosage <- genotypes$gamma[[2]] + 2 * genotypes$gamma[[3]]
  AF <- rowMeans(dosage) / 2
  AF[is.na(AF)] <- 0
  R2 <- apply(dosage, 1, var) / (2 * AF * (1 - AF))
  R2[is.na(R2)] <- 0
  R2_constrained <- sapply(R2, constrain_R2)
  hwe_p <- getHWE(dosage2hard(genotypes))
  filter_AF <- sapply(AF, function(x) filter_by_AF(x, MAF_cutoff))
  filter_R2 <- sapply(R2, function(x) filter_by_R2(x, R2_cutoff_up, R2_cutoff_down))
  filter_HWE <- sapply(hwe_p, function(x) filter_by_HWE(x, 1e-6))
  filter <- paste(filter_AF, filter_R2, filter_HWE, sep=",")
  filter[filter=="PASS,PASS,pASS"] <- "PASS"
  filter[filter!="PASS,PASS,pASS"] <- gsub(",PASS", "", gsub("PASS,", "", filter[filter!="PASS,PASS,pASS"]))

  ## Write into a VCF file
  if(vcf){
    samplelist <- paste(colnames(dosage), collapse="\t")
    header <- paste(
      "##fileformat=VCFv4.2",
      "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
      "##INFO=<ID=R2,Number=1,Type=Float,Description=\"R-square, encoded as var(G)/2p(1-p), where G is dosage genotype and p is allele frequency. Variants with 1<R2<=1.1 are constrained to 1. Variants with R2>1.1 (marked as .) are recommended to remove.\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Genotype dosage.\">",
      "##FORMAT=<ID=RAI,Number=1,Type=Float,Description=\"RAI (Ratio of Alternative allele Intensity).\">",
      "##FORMAT=<ID=GP,Number=1,Type=Float,Description=\"Genotype probability of reference homozygous, heterozygous, and alternative homozygous produced by using ewastools::call_genotypes.\">",
      paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", samplelist),
      sep = "\n"
    )
    hardgeno <- dosage2hard(genotypes)
    geno <- matrix(nrow=nrow(dosage), ncol=ncol(dosage))
    colnames(geno) <- colnames(dosage)
    for(i in 1:nrow(geno)){
      for(j in 1:ncol(geno)){
        geno[i,j] <- paste0(
          hardgeno[i,j], ":",
          round(dosage[i,j], 2), ":", 
          round(genotypes$snps[i,j], 2), ":", 
          round(genotypes$gamma[[1]][i,j], 2), ",", 
          round(genotypes$gamma[[2]][i,j], 2), ",", 
          round(genotypes$gamma[[3]][i,j], 2))
      }
    }
    geno <- as.data.frame(cbind(CpG = rownames(dosage), geno))
    if(type=="snp_probe"){
      data(probeInfo_snp); probeInfo = probeInfo_snp
    }else if(type=="ccs_snp_probe"){
      data(probeInfo_ccs); probeInfo = probeInfo_ccs
    }else if(type=="typeII_probe"){
      data(probeInfo_typeII); probeInfo = probeInfo_typeII
    }else{
      print("Error: misspecified probe types!")
      return
    }
    vcf <- cbind(
      probeInfo[probeInfo$CpG %in% rownames(genotypes$snps), 1:6], 
      tibble(QUAL=".", FILTER=filter, INFO=paste0("AF=", round(AF, 3), ";R2=", R2_constrained), FORMAT="GT:DS:RAI:GP")
    ) %>% left_join(geno, by=c("CpG"))
    vcf <- vcf[, -6]
    write.table(header, file=paste0("genotypes.", type, ".vcf"), sep="\t", row.names=F, quote=F, col.names=F)
    write.table(vcf, file=paste0("genotypes.", type, ".vcf"), sep="\t", row.names=F, quote=F, col.names=F, append=T)
  }
  ## Filter dosage
  dosage <- dosage[filter=="PASS",,drop=F]
  dosage
}
