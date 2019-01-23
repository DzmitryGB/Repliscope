#' Sequence read coverage for wild type replicating sample
#'
#' Sequence read coverage for wild type replicating sample 
#' (T7107 strain). The cells were stained with DNA dye and sorted
#' based on DNA content into S phase fraction. Extracted DNA 
#' was sequenced and mapped to sacCer3 genome. Unique reads were
#' calculated in 1 kb genomic bins.
#'
#' @docType data
#'
#' @usage data(W303_S)
#'
#' @format data frame with 11820 rows and 5 variables:
#' \describe{
#'		\item{chrom}{short chromosome name}
#'		\item{chromStart}{left chromosome coordinate}
#'		\item{chromEnd}{right chromosome coordinate}
#'		\item{name}{sample name}
#'		\item{score}{read number in current bin}
#'	}
#'
#' @keywords datasets sortSeq bed coverage
#'
#' @references Natsume et al. (2013) Mol Cell 50(5):661-74
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23746350}{PubMed})
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/sra/SRX204358[accn]}{SRA}
#'
#' @examples
#' data(W303_S)
"W303_S"
