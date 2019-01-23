#' Sequence read coverage for wild type S.cerevisiae W303 non-replicating sample.
#'
#' Sequence read coverage for wild type non-replicating sample 
#' (T7107 strain). The cells were stained with DNA dye and sorted
#' based on DNA content into G2/M phase fraction. Extracted DNA 
#' was sequenced and mapped to sacCer3 genome. Unique reads were
#' calculated in 1 kb genomic bins.
#'
#' @docType data
#'
#' @usage data(W303_G2)
#'
#' @format data frame with 11350 rows and 5 variables:
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
#' @source \href{https://www.ncbi.nlm.nih.gov/sra/SRX204357[accn]}{SRA}
#'
#' @examples
#' data(W303_G2)
"W303_G2"
