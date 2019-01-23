#' A function to remove outliers from the "ratio" column of a supplied ratio dataframe
#' trimRatio is applied to the calculated ratio of read counts from a replicating to
#' a non-replicating samples.
#' @param ratioDF A ratio dataframe containing 'ratio' column (dataframe, required).
#' @param loLim Low limit threshold (double, required).
#' @param hiLim High limit threshold (double, required).
#' @keywords replication genomics bioinformatics
#' @export
#' @examples
#' W303 <- trimRatio(W303,0.5,1.5)

trimRatio <- function(ratioDF,loLim,hiLim) {
  loLim <- as.numeric(loLim)
  hiLim <- as.numeric(hiLim)
  if (hiLim > loLim) {
    if ("tmpRatio" %in% colnames(ratioDF)) {
      ratioDF <- ratioDF[!(ratioDF$tmpRatio<loLim | ratioDF$tmpRatio>hiLim),]
    } else {
      ratioDF <- ratioDF[!(ratioDF$ratio<loLim | ratioDF$ratio>hiLim),]
    }
  } else if (hiLim == loLim) {
    stop("Supplied limits are equal. I quit!")
  } else {
    warning("'loLim' is higher than 'hiLim'. Swapping the parameters.")
    if ("tmpRatio" %in% colnames(ratioDF)) {
      ratioDF <- ratioDF[!(ratioDF$tmpRatio<hiLim | ratioDF$tmpRatio>loLim),]
    } else {
      ratioDF <- ratioDF[!(ratioDF$ratio<hiLim | ratioDF$ratio>loLim),]
    }
  }
  return(ratioDF)
}
