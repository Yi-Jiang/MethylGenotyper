% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{probeInfo_snp}
\alias{probeInfo_snp}
\title{SNP probe information}
\format{
A data frame with 53 rows and 8 columns:
\describe{
\item{Chr}{Chromosome ID}
\item{Pos}{Position}
\item{SNP}{SNP ID targeted by the CpG}
\item{RefAllele}{Reference allele}
\item{AltAllele}{Alternative allele}
\item{CpG}{CpG}
\item{Color}{Color channel}
\item{Group}{Probe types, color channel, and signal corresponds to alternative allele}
}
}
\source{
\url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
}
\usage{
data(probeInfo_snp)
}
\description{
A dataset containing SNP probe information. Only autosome probes are included.
}
\keyword{datasets}
