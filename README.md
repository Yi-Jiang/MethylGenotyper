# MethylGenotyper: Accurate estimation of SNP genotypes from DNA methylation data

## Description
`MethylGenotyper` aims to call genotypes from three types of probes in Infinium methylation array: SNPs probes, Type I probes with SNPs at the extension bases causing color-channel-switching, and Type II probes with SNPs at the extension bases. By properly modeling the relationship between methylation intensity and SNP genotypes, **MethylGenotyper can produce accurate genotypes for ~4000 SNPs from EPIC array and ~2000 SNPs from 450K array, enables accurate estimation of population structure and genetic relatedness**. We have wrapped all these functions into MethylGenotyper and strongly encourage researchers to incorporate it into the EWAS pipeline to maximize statistical power and avoid spurious association signals.

## Vignette and reference manual
Please read the vignette for detailed description and recommended workflow at [vignette/vignette.pdf](vignette/vignette.pdf) and the reference manual at [vignette/MethylGenotyper_1.0.pdf](vignette/MethylGenotyper_1.0.pdf).

## Installation
```R
library(devtools)
install_github("Yi-Jiang/MethylGenotyper")
```
