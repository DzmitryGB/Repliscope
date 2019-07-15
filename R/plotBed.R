#' A function to boxplot 'score' column of a BED dataframe,
#' per unique chromosome name in the 'chrom' column. The resulting plot also highlights outliers
#' based on the inter quartile range (IQR). The genome wide median is plotted as a pink line through the boxplots.
#' @param bed A dataframe containing 'score' and 'chrom' columns (dataframe).
#' @param plotting Should the plot object be sent to the default device? (boolean, defaults to TRUE).
#' @keywords boxplot BED genomics bioinformatics
#' @import ggplot2
#' @export
#' @examples
#' plotBed(W303_S)
#' plotObject <- plotBed(W303_G2,plotting=FALSE)

plotBed <- function(bed,plotting=TRUE) {
  ## the line below is to appease CRAN checks as it doesn't understand data.frame's variables
  chrom <- score <- NULL
  ## plotting theme
  theTheme<-theme(
    #text = element_text(family="Arial"),
    plot.title = element_text(size=rel(1.6),hjust=0.5,face="bold"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.75),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.5)),
    axis.text.x = element_text(angle=45,hjust=1),
    legend.position='none',
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
  )
  nChr <- length(unique(bed$chrom))
  scoreStats <- summary(bed$score)
  maxScore <- scoreStats[6]
  minScore <- scoreStats[1]
  iqr <- as.numeric(scoreStats[5]-scoreStats[2])
  top3 <- as.numeric(scoreStats[5]+3*iqr)
  low3 <- scoreStats[2]-3*iqr
  boxplot <- ggplot(bed,aes(chrom,score))
  if (as.numeric(scoreStats[6]) > as.numeric(scoreStats[5]+3*iqr)) {
    boxplot <- boxplot + geom_rect(
      xmin=1-nChr*0.03,
      ymin=top3,
      xmax=nChr+nChr*0.03,
      ymax=maxScore+(maxScore-minScore)*0.02,
      fill="gray95",
      size=0
    )
  }
  if (minScore < low3) {
    boxplot <- boxplot + geom_rect(
      xmin=1-nChr*0.03,
      ymin=minScore-(maxScore-minScore)*0.02,
      xmax=nChr+nChr*0.03,
      ymax=low3,
      fill="gray95",
      size=0
    )
  }
  if (bed$name[1]=="name" | bed$name[1]=="") {
    title <- "Read count distribution"
  } else {
    title <- paste0("Read count distribution (",bed$name[1],")")
  }
  boxplot <- boxplot + geom_boxplot()
  boxplot <- boxplot + geom_hline(yintercept=as.numeric(scoreStats[3]),size=.75,color="orchid")  ##  Median line
  boxplot <- boxplot + coord_cartesian(xlim=c(1,nChr),ylim=c(minScore,maxScore))
  boxplot <- boxplot + ggtitle(title)
  boxplot <- boxplot + labs(y="Score value")
  boxplot <- boxplot + theTheme
  if (plotting) { boxplot } else { return(boxplot) }
}
