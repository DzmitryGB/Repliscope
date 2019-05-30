#' A function to normalise ratio values from 'ratio' column of the provided dataframe
#' to fit biologically-relevant scale. It scales values
#' either using supplied 'rFactor' value or automatically to best fit 1 to 2 scale
#' (the upper limit of the scale may be adjusted with the upperLimit parameter).
#' Normalisation factor used is stored in 'ratioFactor' column and also passed as the
#' dataframe comment. To extract it, use 'attributes(mergedBed)$comment'
#' @aliases normalizeRatio
#' @param ratioDF A ratio dataframe containing 'ratio' column (dataframe).
#' @param rFactor Value to normalise by, related to replication progression (numeric, optional).
#' @param upperLimit Top value for the scale, defaults to 2 (numeric, optional, defaults).
#' @param replace Should the existing 'ratio' values be overwritten or stored in a new column
#'  (boolean, defaults to TRUE).
#' @export
#' @importFrom stats optimise
#' @keywords replication genomics bioinformatics
#' @examples
#' ratioDF <- normaliseRatio(W303) ## scales to 1 to 2 range, replaces original values.
#' ratioDF <- normaliseRatio(W303,rFactor=1.41,replace=FALSE)
#' # (multiplies score values by 1.41 and keeps the original values)

normaliseRatio <- function(ratioDF,rFactor=NULL,upperLimit=2,replace=TRUE) {
  if (is.null(rFactor)) {
    if (upperLimit > 1) {
      f <- function (factor,values,ceiling=upperLimit) abs(sum(factor*values[factor*values>ceiling]-ceiling)-sum(1-factor*values[factor*values<1]))
      xmin <- optimise(f,c(1,upperLimit),tol = 0.001,values=na.omit(ratioDF$ratio),maximum=F)
      ratioFactor <- xmin[[1]]
    } else { stop("The upperLimit must be bigger than 1!") }
  } else {
    ratioFactor <- as.numeric(rFactor)
  }
  if (replace) {
    ratioDF$ratio <- round(ratioDF$ratio*ratioFactor,3)
  } else {
    ratioDF$tmpRatio <- ratioDF$ratio*ratioFactor
  }
  ratioDF$ratioFactor <- round(ratioFactor*as.numeric(ratioDF$ratioFactor),3)
  comment(ratioDF) <- as.character(round(ratioFactor,3))
  return(ratioDF)
}
