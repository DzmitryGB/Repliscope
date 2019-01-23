#' S.cerevisiae genome information  
#'
#' sacCer3 meta information: chromosome sizes, 
#' centromere and replication origin positions.
#'
#' @encoding UTF-8
#'
#' @docType data
#'
#' @usage data(sacCer3)
#'
#' @format List containing three dataframes
#'
#'		\describe{
#'			\item{genome}{Chromosome information dataframe}
#'			\item{cen}{Centromere information dataframe}
#'			\item{ori}{Replication origin information dataframe}
#'		}
#'
#' @keywords datasets sacCer3 replication
#'
#' @references Siow et al. (2011) NAR 40(Database issue):D682-6
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/22121216}{PubMed})
#'
#' @source Replication origin information: (\href{http://cerevisiae.oridb.org/}{OriDB})
#'
#' @examples
#' data(sacCer3)
"sacCer3"
