#' A function to compare two replication profiles
#'
#' compareRatios takes two ratio dataframes that were binned the same way and uses z-score
#' statistics to find p-values of their differences. The function outputs a combined dataframe containing
#' the two input ratio dataframes in a long format with added 'p.value' column.
#' @param ratio1 Ratio dataframe, or a string containing name of a ratio dataframe (dataframe or string).
#' @param ratio2 Ratio dataframe, or a string containing name of a ratio dataframe (dataframe or string).
#' @keywords BED genomics bioinformatics
#' @export
#' @importFrom stats lm na.omit pnorm quantile qqnorm
#' @examples
#' ratioDFs <- compareRatios(W303norm,Dbf4myc)

compareRatios <- function(ratio1,ratio2) {
  if (typeof(ratio1)=="character" & typeof(ratio2)=="character") {
    nameRatio1<-ratio1
    nameRatio2<-ratio2
    try(ratio1 <- eval(as.name(nameRatio1)))
    try(ratio2 <- eval(as.name(nameRatio2)))
  } else {
    nameRatio1<-deparse(substitute(ratio1))
    nameRatio2<-deparse(substitute(ratio2))
  }
  ## the line below is to appease CRAN checks as it doesn't understand data.frame's variables
  name.rep.x <- name.rep.y <- chrom <- NULL
  ## calculations
  if (typeof(ratio1) != "list") { stop("Could not find dataframe ",nameRatio1,"!")
  } else if (typeof(ratio2) != "list") { stop("Could not find dataframe ",nameRatio2,"!")
  } else {
    ratios <- merge(ratio1,ratio2,by=c("chrom","chromStart","chromEnd"),all=T,sort=F)
    #~ 		cat("\tComparing replication profiles between ",nameRatio1," and ",nameRatio2,"\n",sep="")
    calc <- data.frame(diff=ratios$ratio.x-ratios$ratio.y,means=(ratios$ratio.x+ratios$ratio.y)/2)
    q=qqnorm(calc$diff, plot.it=F)
    qnt_x <- q$x[! quantile(q$x,0.25,na.rm=T)>=q$x & q$x<=quantile(q$x,0.75,na.rm=T)]
    qnt_y <- q$y[! quantile(q$x,0.25,na.rm=T)>=q$x & q$x<=quantile(q$x,0.75,na.rm=T)]
    my_fit=lm(qnt_y ~ qnt_x)
    est_mean=my_fit$coefficients[1]
    est_sd=my_fit$coefficients[2]
    calc$z.score <- (abs(calc$diff-est_mean))/est_sd
    ratios$p.value <- 1-pnorm(calc$z.score)
    # subset DF by x and y
    longRatio.x <- subset(ratios,!is.na(name.rep.x),select=names(ratios)[c(1:3,grep(".x$",names(ratios)),ncol(ratios))])
    longRatio.y <- subset(ratios,!is.na(name.rep.y),select=names(ratios)[c(1:3,grep(".y$",names(ratios)),ncol(ratios))])
    # rename subsetted DFs
    names(longRatio.x) <- gsub("\\.x$","",names(longRatio.x))
    names(longRatio.y) <- gsub("\\.y$","",names(longRatio.y))
    longRatio.y$p.value <- NA
    # reorder y DF
    chroms <- unique(ratio2$chrom)
    orderedDF <- longRatio.y[0,]
    for (chr in chroms) {
      tmp <- subset(longRatio.y,chrom==chr)
      orderedDF <- rbind(orderedDF,tmp[order(tmp$chromStart),])
    }
    orderedDF$chrom <- factor(orderedDF$chrom, levels=chroms)
    longRatios <- rbind(longRatio.x,orderedDF)
    rownames(longRatios) <- 1:nrow(longRatios)
    longRatios$ratioFactor <- factor(longRatios$ratioFactor,levels=unique(longRatios$ratioFactor))
    return(longRatios)
  }
}
