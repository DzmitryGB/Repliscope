#' A function to calculate Trep values from a sync-seq experiment
#' calcTrep function fits a Boltzman sigmoid function into relative copy number datapoints for every
#' genomic bin of the provided sync-seq merged dataframe. It then extracts time at which half of the
#' cells have this genomic bin replicated (Trep). The output of the function is a dataframe containing
#' Trep and TrepErr data for every genomic bin in a BED-like format.
#'
#' @param ratioDFs A merged ratios dataframe containing sync-seq samples (dataframe).
#' @param times Time series data in the same order as in the ratioDFs (numeric vector).
#' @export
#' @importFrom stats nls
#' @keywords trep sigmoid replication genomics bioinformatics
#' @examples
#' TrepDF <- calcTrep(subset(syncSeq[["data"]],chrom=="chrI"),times=c(25,30,35,40,45,50,90))

calcTrep <- function(ratioDFs,times) {
  # get time series limits to do sanity check later
  minTime <- min(times)
  maxTime <- max(times)
  # add 0 timepoint to the supplied time series
  x<-append(0,times)
  # supply constants for the Boltzman function
  a <- 2 # max value of y
  b <- 1 # min value of y
  # find different samples in the supplied ratioDFs:
  samples <- unique(ratioDFs[,c("name.rep","name.nonRep")])
  rownames(samples) <- 1:nrow(samples)
  i <- dim(samples)[1]
  if (i < 2) stop("More than one sample is required for Trep calculation. Use 'rbind()' to merge multiple samples.")
  # get values based on unique coordinates
  resultDF <- data.frame("chrom"=character(),"chromStart"=character(),"chromEnd"=character(),"Trep"=character(),"TrepErr"=character(),stringsAsFactors=F)
  binsDF <- unique(ratioDFs[,c("chrom","chromStart","chromEnd")])
  for (binRow in 1:nrow(binsDF)) {
    bin <- binsDF[binRow,]
    ratioValues <- numeric() ## initiate empty vector to store all ratios for a bin
    # extract available ratios, filling with NA in place of missing values
    for (j in 1:i) {
      aSample <- samples[j,]
      ratioValueRow <- which(ratioDFs$name.rep==aSample[1,1] & ratioDFs$name.nonRep==aSample[1,2] & ratioDFs$chrom==bin[1,1] & ratioDFs$chromStart==bin[1,2] & ratioDFs$chromEnd==bin[1,3])
      if (as.logical(length(ratioValueRow[1] == 1))) {
        ratioValues <- append(ratioValues,ratioDFs$ratio[ratioValueRow])
      } else if (as.logical(length(ratioValueRow[1] == 0))) {
        ratioValues <- append(ratioValues,NA)
      } else {stop("There are non-unique samples in the supplied dataframe. I quit!")}
    }
    y<-append(1,ratioValues)
    # fit sigmoidal function to the extracted values, where c is Trep and d is slope
    boltzmann <- NULL
    try(boltzmann <- nls(y ~ a+(b-a)/(1+exp((x-c)/d)), start=list(c=25,d=1)),silent=T)
    if (!is.null(boltzmann)) {
      params <- summary(boltzmann)
      Trep <- as.character(round(params$parameters[1,1],1))
      TrepErr <- as.character(round(params$parameters[1,2],2))
    } else {
      Trep <- TrepErr <- NA
    }
    # write into results DF
    resultDF[nrow(resultDF)+1,] <- c(as.character(bin[1,1]),as.character(bin[1,2]),as.character(bin[1,3]),Trep,TrepErr)
  }
  resultDF$chrom <- factor(resultDF$chrom,levels=unique(resultDF$chrom))
  resultDF$chromStart <- as.integer(resultDF$chromStart)
  resultDF$chromEnd <- as.integer(resultDF$chromEnd)
  resultDF$Trep <- as.numeric(resultDF$Trep)
  resultDF$TrepErr <- as.numeric(resultDF$TrepErr)
  # sanity check - the Trep values must be within the time series limits
  bogusRows <- which(resultDF$Trep < minTime | resultDF$Trep > maxTime)
  if (length(bogusRows)[1] != 0) {
    resultDF[bogusRows,"Trep"] <- NA
    resultDF[bogusRows,"TrepErr"] <- NA
  }
  return(resultDF)
}
