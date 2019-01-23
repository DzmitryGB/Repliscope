#' A helper function to create axis ticks and human readable labels.
#'
#' makeLabels is called by plotGenome() and plotCoverage() functions. It creates a dataframe
#' containing two columns: 'ticks' and 'labels'. 'Ticks' contains axis ticks coordinates, 'labels'
#' will contain human readable lables for the ticks (using prefixes and optional units).
#' @param theMin Minimum value for the scale (double).
#' @param theMax Maximum value for the scale (double).
#' @param unit Unit to use for the labels (string, optional)
#' @keywords plotting genomics bioinformatics
#' @export
#' @examples
#' labels <- makeLabels(0,1200000,"b")

makeLabels <- function(theMin,theMax,unit="") {
  if (!is.na(as.numeric(theMin)) & !is.na(as.numeric(theMax))) {
    step <- 10^floor(log10(theMax-theMin))
    minorFactor <- 2
    if ((theMax-theMin)/step < 6) {
      step <- step/2
      minorFactor <- 2
    }
    if ((theMax-theMin)/step < 6) {
      step <- step/2
      minorFactor <- 2.5
    }
    if ((theMax-theMin)/step < 6) {
      step <- step/2.5
      minorFactor <- 2
    }
    minTick <- step*ceiling(theMin/step)
    if ((minTick-theMin)/step > 0.5) minTick <- minTick-step
    maxTick <- step*floor(theMax/step)
    if ((theMax-maxTick)/step > 0.5) maxTick <- maxTick+step
    labels <- seq(from = minTick, to = maxTick, by = step)
    if (log10(step) < 0) {
      prefix <- switch(abs(ceiling(log10(step)/3))+1,"","m","\u00B5","n","p","f","a")
      divisor <- switch(abs(ceiling(log10(step)/3))+1,1,0.001,0.000001,0.000000001,0.000000000001,0.000000000000001,0.000000000000000001)
    } else {
      prefix <- switch(floor(log10(step)/3)+1,"","k","M","G","T","P","E")
      divisor <- switch(floor(log10(step)/3)+1,1,1000,1000000,1000000000,1000000000000,1000000000000000,1000000000000000000)
    }
    names <- paste0(labels/divisor,prefix,unit)
    names <- gsub(paste0("^",0,prefix,unit),"0",names)
    labels <- data.frame(ticks=labels,labels=names)
    attr(labels,"minorTickFactor") <- as.character(minorFactor)
    return(labels)
  } else {
    stop("There is a problem with the supplied minimum and maximum!")
  }
}
