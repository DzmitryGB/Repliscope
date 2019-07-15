#' Trep data calculated from syncSeq[["data"]]
#'
#' Trep is median replication time, expressed in minutes after release G1 arrest.
#' It is calculated from multiple relative copy number datapoints across timeseries
#' of a cell cycle experiment. For every genomic bin, a sigmoid function is fitted
#' and its midpoint is reported.
#'
#' @docType data
#'
#' @usage data(TrepDF)
#'
#' @format data frame with 11341 rows and 5 variables:
#' \describe{
#'		\item{chrom}{short chromosome name}
#'		\item{chromStart}{left chromosome coordinate}
#'		\item{chromEnd}{right chromosome coordinate}
#'		\item{Trep}{calculated Trep value}
#'		\item{TrepErr}{error from sigmoid function fitting}
#'	}
#'
#' @keywords datasets syncSeq Trep
#'
#' @references Müller et al. (2014) NAR 42(1):e3
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24089142}{PubMed})
#'
#' @examples
#' data(TrepDF)
"TrepDF"
