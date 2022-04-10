
#' EPIC manifest file
#' 
#' A dataset containing all EPIC probes information.
#' 
#' @usage data(mnfst)
#' @format A data frame with 866554 rows and 5 columns:
#' \describe{
#'  \item{Name}
#'  \item{AddressA_ID}
#'  \item{AddressB_ID}
#'  \item{Infinium_Design_Type}
#'  \item{Color_Channel}
#' }
#' @source \url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
"mnfst"

#' Genotyping probe information
#' 
#' A dataset containing genotyping probe information.
#' 
#' @usage data(probeInfo_geno)
#' @format A data frame with 53 rows and 8 columns:
#' \describe{
#'  \item{Chr}{Chromosome ID}
#'  \item{Pos}{Position}
#'  \item{SNP}{SNP ID targeted by the CpG}
#'  \item{RefAllele}{Reference allele}
#'  \item{AltAllele}{Alternative allele}
#'  \item{CpG}{CpG}
#'  \item{Color}{Color channel}
#'  \item{Group}{Probe types, color channel, and signal corresponds to alternative allele}
#' }
#' @source \url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
"probeInfo_geno"

#' CCS SNP probe information
#' 
#' A dataset containing CCS SNP probe information.
#' 
#' @usage data(probeInfo_ccs)
#' @format A data frame with 132 rows and 8 columns:
#' \describe{
#'  \item{Chr}{Chromosome ID}
#'  \item{Pos}{Position}
#'  \item{SNP}{SNP ID targeted by the CpG}
#'  \item{RefAllele}{Reference allele}
#'  \item{AltAllele}{Alternative allele}
#'  \item{CpG}{CpG}
#'  \item{Color}{Color channel}
#'  \item{Group}{NA}
#' }
#' @source \url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
"probeInfo_ccs"

#' Type-II SNP probe information
#' 
#' A dataset containing Type-II SNP probe information.
#' 
#' @usage data(probeInfo_ccs)
#' @format A data frame with 4872 rows and 8 columns:
#' \describe{
#'  \item{Chr}{Chromosome ID}
#'  \item{Pos}{Position}
#'  \item{SNP}{SNP ID targeted by the CpG}
#'  \item{RefAllele}{Reference allele}
#'  \item{AltAllele}{Alternative allele}
#'  \item{CpG}{CpG}
#'  \item{Color}{Color channel}
#'  \item{Group}{Probe types, color channel, and signal corresponds to alternative allele}
#' }
#' @source \url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
"probeInfo_typeII"

#' Probe list
#' 
#' A dataset containing the list of 53 genotyping probes on autosomes, 132 CCS SNP probes, and 4872 type-II SNP probes.
#' 
#' @usage data(probelist)
#' @format A data frame with 5057 rows and 2 columns:
#' \describe{
#'  \item{CpG}{CpG list}
#'  \item{Type}{Probe types}
#' }
#' @seealso Probe informations: \code{\link{probeInfo_geno()}}, \code{\link{probeInfo_ccs()}}, and \code{\link{probeInfo_typeII()}}.
"probelist"

