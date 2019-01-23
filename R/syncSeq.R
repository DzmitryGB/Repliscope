#' Replication profiles budding yeast arrest-release samples
#'
#' Replication profiles of wild type S. cerevisiae arrest-release samples 
#' (AUY077 strain). The cells were arrested in G1 with alpha-factor
#' followed by release using pronase. The samples were collected
#' before the release (aFactor) and various time intervals after
#' the release (25min,30min,35min,40min,45min,50min and 90min)
#' Extracted DNA was sequenced and mapped to sacCer3 genome. Unique
#' reads for replicating (post-release) and non-replicating (aFactor)
#' samples were calculated in 1 kb genomic bins. The ratios were 
#' created by dividing 'score' values from replicating samples by 
#' non-replicating sample 'score' values, adjusted by total number
#' of reads. The ratio values were further adjusted based on bulk 
#' genome replication (as determined by flow cytometry), to put the 
#' values onto biologically relevant relative copy number scale from 
#' 1 to 2. The relative copy number values were smoothed using cubic 
#' spline.
#'
#' @docType data
#'
#' @usage data(syncSeq)
#'
#' @encoding UTF-8
#'
#' @format List containing two data frames
#'	\describe{
#'		\item{data}{syncSeq replication profiles data. Columns:
#'			\strong{chrom} (short chromosome name),
#'			\strong{chromStart} (left chromosome coordinate),
#'			\strong{chromEnd} (right chromosome coordinate),
#'			\strong{name.rep} (replicating sample name),
#'			\strong{name.nonRep} (non-replicating sample name),
#'			\strong{ratio} (ratio value in the current bin),
#'			\strong{ratioFactor} (adjustment factor used for the current ratio),
#'			\strong{group} (Group number of the current bin),
#'			\strong{splineSmooth} (Smoothed ratio value)
#'		}
#'		\item{guide}{Guide dataframe for plotting the syncSeq data
#'			\strong{order} (Order to plot data in),
#'			\strong{name.rep} (Name of replicating sample),
#'			\strong{name.nonRep} (Name of non-replicating sample),
#'			\strong{raw} (Should raw data be plotted?),
#'			\strong{smooth} (Should smooth data be plotted?),
#'			\strong{color} (Color to plot the profile in)
#'		}
#'	}
#'
#' @keywords datasets syncSeq replication
#'
#' @references Müller et al. (2014) NAR 42(1):e3
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/24089142}{PubMed})
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48212}{GEO}
#'
#' @examples
#' data(syncSeq)
"syncSeq"
