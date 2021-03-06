% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccs_snp_probe.R
\name{getRAI_ccs}
\alias{getRAI_ccs}
\title{Get RAI (Ratio of Alternative allele Intensity) for CCS probes}
\usage{
getRAI_ccs(rgData, pop = "ALL")
}
\arguments{
\item{rgData}{Noob and dye-bias corrected signals produced by using \code{correct_noob_dye}.}

\item{pop}{Population. One of EAS, AMR, AFR, EUR, SAS, and ALL. Only probes with MAF of matching population > 0.01 will be kept. Only effective when train=TRUE.}
}
\value{
RAI (Ratio of Alternative allele Intensity).
}
\description{
Get RAI (Ratio of Alternative allele Intensity) for CCS probes
}
