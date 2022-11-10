
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

#' SNP probe information
#' 
#' A dataset containing SNP probe information. Only autosome probes are included.
#' 
#' @usage data(probeInfo_snp)
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
"probeInfo_snp"

#' Type I CCS probe information
#' 
#' A dataset containing Type I CCS probe information.
#' 
#' @usage data(probeInfo_typeI)
#' @format A data frame with 717 rows and 16 columns:
#' \describe{
#'  \item{Chr}{Chromosome ID}
#'  \item{Pos}{Position}
#'  \item{SNP}{SNP ID targeted by the CpG}
#'  \item{RefAllele}{Reference allele}
#'  \item{AltAllele}{Alternative allele}
#'  \item{CpG}{CpG}
#'  \item{Color}{Color channel}
#'  \item{Group}{NA}
#'  \item{ALL_AF}{Allele frequency of all population}
#'  \item{EAS_AF}{Allele frequency of East Asian}
#'  \item{AMR_AF}{Allele frequency of American}
#'  \item{AFR_AF}{Allele frequency of African}
#'  \item{EUR_AF}{Allele frequency of European}
#'  \item{SAS_AF}{Allele frequency of South Asian}
#'  \item{h_0.1}{Passed peak density test or not}
#'  \item{loc_pass}{Passed peak position test or not}
#' }
#' @source \url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
"probeInfo_typeI"

#' Type II CCS probe information
#' 
#' A dataset containing information of Type II CCS probes with SNPs at the extension bases. We only consider the situation that the alternative allele is A/T and the reference allele is C/G.
#' 
#' @usage data(probeInfo_typeII)
#' @format A data frame with 26420 rows and 16 columns:
#' \describe{
#'  \item{Chr}{Chromosome ID}
#'  \item{Pos}{Position}
#'  \item{SNP}{SNP ID targeted by the CpG}
#'  \item{RefAllele}{Reference allele}
#'  \item{AltAllele}{Alternative allele}
#'  \item{CpG}{CpG}
#'  \item{Color}{Color channel}
#'  \item{Group}{Probe types, color channel, and signal corresponds to alternative allele}
#'  \item{ALL_AF}{Allele frequency of all population}
#'  \item{EAS_AF}{Allele frequency of East Asian}
#'  \item{AMR_AF}{Allele frequency of American}
#'  \item{AFR_AF}{Allele frequency of African}
#'  \item{EUR_AF}{Allele frequency of European}
#'  \item{SAS_AF}{Allele frequency of South Asian}
#'  \item{h_0.1}{Passed peak density test or not}
#'  \item{loc_pass}{Passed peak position test or not}
#' }
#' @source \url{https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip}
"probeInfo_typeII"

#' Probe list
#' 
#' A dataset containing the list of 53 SNP probes on autosomes, 717 Type I CCS probes, and 27879 type II CCS probes.
#' 
#' @usage data(probelist)
#' @format A data frame with 27190 rows and 2 columns:
#' \describe{
#'  \item{CpG}{CpG list}
#'  \item{Type}{Probe types}
#' }
#' @seealso Probe informations: \code{\link{probeInfo_snp()}}, \code{\link{probeInfo_typeI()}}, and \code{\link{probeInfo_typeII()}}.
"probelist"

