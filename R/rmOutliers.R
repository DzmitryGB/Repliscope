#' A function to remove outliers from the "score" column of a supplied bed dataframe
#' There are three methods: max, IQR and median. Max is used to remove 1 or more maximum
#' values; IQR uses interquartile range to detect outliers, while median method can be
#' used to remove data based on genome-wide median.
#' @param bed A dataframe containing 'score' column (dataframe, required).
#' @param method Method to detect outliers: "max", "IQR" or "median" (string).
#' @param n Number of max values to remove (integer,defaults to 1). Use with "max" method.
#' @param range Number of IQR above the 3rd or below the 1st IQR to set the threshold (double, defaults to 3).
#'  Use with "IQR" method.
#' @param loLim Low limit for the median method (double, defaults to 0.25).
#' @param hiLim High limit for the median method (double).
#' @keywords outliers BED genomics bioinformatics
#' @export
#' @examples
#' bedDF <- rmOutliers(W303_S,method="max",n=2) ## removes 2 rows of data containing 3 top values
#' bedDF <- rmOutliers(W303_S,method="IQR",range=3) ## removes datapoints outside 3 x IQR above the 3rd
#'  # and below the 1st IQR.
#' bedDF <- rmOutliers(W303_S,method="median",loLim=0.25,hiLim=2) # removes datapoints that are lower
#'  # than 0.25 x genome median or above 2 x genome median.

rmOutliers <- function(bed,method,n=1,range=3,loLim=0.25,hiLim=NULL) {
  if (method == "max") {
    if (n!=1) {
      n <- as.integer(n)
    }
    if (n > 0) {
      for (times in 1:n) {
        bed<-bed[-which(bed$score==max(bed$score)),]
      }
    } else {
      stop("'n' must be a positive integer.")
    }
  } else if (method=="IQR") {
    if (range!=3) {
      range <- as.numeric(range)
      if (range <= 0) {
        stop("'range' must be positive number!")
      } else if (range <= 1.5) {
        warning("Those are not likely to be outliers!")
      } else if (range <= 3) {
        warning("Removing suspected outliers!")
      }
    }
    scoreStats <- summary(bed$score)
    iqr <- scoreStats[5]-scoreStats[2]
    bed<-bed[!(bed$score<scoreStats[2]-range*iqr | bed$score>scoreStats[5]+range*iqr),]
  } else if (method=="median") {
    if (!(is.null(loLim) & is.null(hiLim))) {
      median <- as.numeric(summary(bed$score)[3])
      if (!is.null(loLim)) {
        if (!is.na(as.numeric(loLim))) {
          lo <- as.numeric(as.numeric(loLim)*median)
          bed <- bed[bed$score>lo,]
        }
      }
      if (!is.null(hiLim)) {
        if (!is.na(as.numeric(hiLim))) {
          hi <- as.numeric(as.numeric(hiLim)*median)
          bed <- bed[bed$score<hi,]
        }
      }
    } else {
      stop("You must provide either 'loLim' or 'hiLim' for the 'median' method.")
    }
  } else {
    stop("A method must be provided: 'max','IQR' or 'median'!")
  }
  rownames(bed) <- 1:nrow(bed)
  return(bed)
}
