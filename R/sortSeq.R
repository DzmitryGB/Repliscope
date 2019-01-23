#' Replication profiles for wild type and Dbf4-9myc S.cerevisiae samples
#'
#' Replication profiles for wild type and Dbf4-9myc samples 
#' (T7107 and T9394 strains). The cells were stained with DNA dye 
#' and sorted based on DNA content into S or G2/M phase fractions.
#' Extracted DNA was sequenced and mapped to sacCer3 genome. Unique
#' reads for replicating (S) and non-replicating (G2/M) samples were
#' calculated in 1 kb genomic bins. The ratio was created by
#' dividing 'score' values from replicating sample by non-
#' replicating sample 'score' values, adjusted by total number
#' of reads. The ratio values were further adjusted by multiplying
#' them by 1.41 and 1.402 for wild type and Dbf4-9myc samples, 
#' respectively, to put the values onto biologically relevant
#' relative copy number scale from 1 to 2. The relative copy number
#' values were smoothed using cubic spline and compared using z
#' score statistics.
#'
#' @docType data
#'
#' @usage data(sortSeq)
#'
#' @format data frame with 22696 rows and 10 variables:
#' \describe{
#'		\item{chrom}{short chromosome name}
#'		\item{chromStart}{left chromosome coordinate}
#'		\item{chromEnd}{right chromosome coordinate}
#'		\item{name.rep}{replicating sample name}
#'		\item{name.nonRep}{non-replicating sample name}
#'		\item{ratio}{ratio value in the current bin}
#'		\item{ratioFactor}{adjustment factor used for the current ratio}
#'		\item{group}{Group number of the current bin}
#'		\item{splineSmooth}{Smoothed ratio value}
#'		\item{p.value}{Significance of ratio difference between Dbf4myc and W303 samples}
#'	}
#'
#' @keywords datasets sortSeq replication
#'
#' @references Natsume et al. (2013) Mol Cell 50(5):661-74
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23746350}{PubMed})
#'
#' @source Dbf4myc S phase sample: \href{https://www.ncbi.nlm.nih.gov/sra/SRX202404[accn]}{SRA};
#' Dbf4myc G2 sample: \href{https://www.ncbi.nlm.nih.gov/sra/SRX202403[accn]}{SRA};
#' W303 S sample: \href{https://www.ncbi.nlm.nih.gov/sra/SRX204358[accn]}{SRA};
#' W303 G2 sample: \href{https://www.ncbi.nlm.nih.gov/sra/SRX204357[accn]}{SRA}
#'
#' @examples
#' data(sortSeq)
"sortSeq"
