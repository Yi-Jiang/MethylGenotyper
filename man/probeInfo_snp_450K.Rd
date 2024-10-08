% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{probeInfo_snp_450K}
\alias{probeInfo_snp_450K}
\title{SNP probe information for 450K}
\format{
A data frame with 57 rows and 14 columns:
\describe{
\item{Chr}{Chromosome ID}
\item{Pos}{Position}
\item{SNP}{SNP ID targeted by the CpG}
\item{RefAllele}{Reference allele}
\item{AltAllele}{Alternative allele}
\item{CpG}{CpG}
\item{Color}{Color channel}
\item{Group}{Probe types, color channel, and signal corresponds to alternative allele}
\item{ALL_AF}{Allele frequency of all population}
\item{EAS_AF}{Allele frequency of East Asian}
\item{AMR_AF}{Allele frequency of American}
\item{AFR_AF}{Allele frequency of African}
\item{EUR_AF}{Allele frequency of European}
\item{SAS_AF}{Allele frequency of South Asian}
}
}
\source{
\url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
}
\usage{
data(probeInfo_snp_450K)
}
\description{
A dataset containing SNP probe information. Only autosome probes are included.
}
\keyword{datasets}
