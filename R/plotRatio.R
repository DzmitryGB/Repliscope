#' A function to plot a histogram of supplied ratio vector
#' plotRatio plots histogram of values in a supplied vector using ggplot2 and
#' highlights interval between 1 and 2 in green.
#' @param ratio A numeric vector containing raw or smoothed ratio values (vector).
#' @param plotting Should the plot object be sent to the default device? (boolean, defaults to TRUE).
#' @keywords histogram BED genomics bioinformatics
#' @import ggplot2
#' @export
#' @examples
#' plotRatio(W303$ratio)
#' plotObject <- plotRatio(W303$ratio,plotting=FALSE)

plotRatio <- function(ratio,plotting=TRUE) {
  theMax <- max(as.numeric(ratio))
  theMin <-min(as.numeric(ratio))

  if (theMin > 0.99) theMin <- 0.99
  if (theMax < 2.01) theMax <- 2.01

  labels <- makeLabels(theMin,theMax)

  theTheme<-theme(
    #text = element_text(family="Arial"),
    plot.title = element_text(size=rel(2),hjust=0.5,face="bold"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=0.75),
    panel.background = element_blank(),
    axis.line = element_blank(),
    axis.title = element_text(size = rel(1.5)),
    axis.text = element_text(size = rel(1.2)),
    legend.position='none',
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
  )
  hist <- ggplot() + aes(ratio)
  hist <- hist + annotate("rect", xmin=1, xmax=2, ymin=0, ymax=Inf, alpha=0.5, fill="lightgreen")
  hist <- hist + geom_histogram(bins=100,fill='gray90',color='gray30')
  hist <- hist + ggtitle("Ratio distribution")
  hist <- hist + scale_x_continuous(
    name="Ratio values",
    breaks=labels$ticks,
    labels=labels$labels,
    limits=c(theMin,theMax),
    expand=c(0.02,0.02))
  hist <- hist + scale_y_continuous(name="Counts within bin")
  hist <- hist + theTheme
  if (plotting) { hist } else { return(hist) }
}
