#' Replication profile for wild type DS2 H.volcanii
#'
#' Replication profile for H.volcanii wild isolate DS2.
#' Genomic DNA for deep sequencing was isolated from 100 ml culture
#' in stationary phase (A650 > 1, DS2_stat sample) or 1 litre in
#' exponential phase (A650 0.1, DS2_exp sample). Unique
#' reads for the two samples were calculated in 1 kb genomic bins
#' using ASM2568v1 genome assembly. The ratio was created by
#' dividing 'score' values from replicating sample by non-
#' replicating sample 'score' values, adjusted by total number
#' of reads. The ratio values were further adjusted by multiplying
#' them by 1.12 to put the values onto biologically relevant
#' relative copy number scale from 1 to 2.
#'
#' @docType data
#'
#' @usage data(MFAseq)
#'
#' @format data frame with 3887 rows and 7 variables:
#' \describe{
#'		\item{chrom}{short chromosome name}
#'		\item{chromStart}{left chromosome coordinate}
#'		\item{chromEnd}{right chromosome coordinate}
#'		\item{name.rep}{replicating sample name}
#'		\item{name.nonRep}{non-replicating sample name}
#'		\item{ratio}{ratio value in the current bin}
#'		\item{ratioFactor}{adjustment factor used for the current ratio}
#'	}
#'
#' @keywords datasets sortSeq replication
#'
#' @references Hawkins et al. (2013) Nature 503(7477):544-547
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24185008}{PubMed})
#'
#' @source DS2_exp exponential phase sample: \href{https://www.ncbi.nlm.nih.gov/sra/SRX202169[accn]}{SRA};
#' DS2_stat stationary sample: \href{https://www.ncbi.nlm.nih.gov/sra/SRX202170[accn]}{SRA}
#'
#' @examples
#' data(MFAseq)
"MFAseq"
