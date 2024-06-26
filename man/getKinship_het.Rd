% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post_genotype_calling.R
\name{getKinship_het}
\alias{getKinship_het}
\title{Get kinship coefficients using the SEEKIN-het estimator}
\usage{
getKinship_het(dosage, indAF)
}
\arguments{
\item{dosage}{A matrix of genotype calls. Provide probes as rows and samples as columns.}

\item{indAF}{A matrix of individual-specific AFs. Provide probes as rows and samples as columns.}
}
\value{
A data frame containing kinship coefficient (Phi) and sample relationships between each two samples.
}
\description{
Only SNPs with missing rate < 10\% were used.
}
