#' A function to scatterplot 'score' column of a BED dataframe
#' plotCoverage function plots values in the ‘score’ column of the supplied bed dataframe as a function of
#' chromosome coordinates. The genome wide median is plotted as a pink line.
#' @param bed A dataframe containing 'score','chrom','chromStart' and 'chromEnd' columns (dataframe).
#' @param region Only plot for the provided region in the format 'chrI:1000-3000' (string, optional).
#' @param plotting Should the plot object be sent to the default device? (boolean, defaults to TRUE).
#' @keywords scatterplot BED genomics bioinformatics
#' @import ggplot2
#' @export
#' @examples
#' plotCoverage(W303_G2)
#' plotObject <- plotCoverage(W303_S,plotting=FALSE)

plotCoverage <- function(bed,region=FALSE,plotting=TRUE) {
  ## the line below is to appease CRAN checks as it doesn't understand data.frame's variables
  chromStart <- chromEnd <- midY <- score <- NULL
  ##  load ggplot2 library and set plotting theme
  fontSize <- 12
  pointSize <- 1 / 0.352777778
  gplotTheme<-theme(
    text = element_text(size=fontSize),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5,size=rel(1.5)),
    axis.line = element_line(colour="gray50",size=.5),
    axis.title = element_text(size = rel(1.5),face="bold"),
    axis.text.x = element_text(size=rel(1.4)),
    legend.position="none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )
  theMedian <- as.numeric(summary(bed$score)["Median"])
  ylims <- c(min(bed$score),max(bed$score))
  ##  Create genome dataframe using a helper function
  if (region!=F) {
    plottingRegion <- T
    region <- as.character(region)
    genome <- makeGenome(bed,region)
    if (length(grep(":",region))==1) {
      region <- unlist(strsplit(region,":"))
      tmp <- unlist(strsplit(region[2],"-"))
      region <- c(region[1],tmp)
      bed<-bed[bed$chrom==region[1] & bed$chromStart>=as.integer(region[2]) & bed$chromEnd<=as.integer(region[3]),]
      #bed<-bed[bed$chromStart>=as.integer(region[2]),]
      #bed<-bed[bed$chromEnd<=as.integer(region[3]),]
    } else {
      bed<-bed[bed$chrom==region]
    }
  } else {
    genome <- makeGenome(bed)
    plottingRegion <- F
  }
  xAxisName <- attributes(genome)$axisName
  ##  Make nice x-axis labels using a helper function
  theMax <- max(genome$chromEnd)
  theMin <- min(genome$chromStart)
  labels <- makeLabels(theMin,theMax,"b")
  #~ 	minorTickFactor <- as.numeric(attributes(labels)$minorTickFactor)
  ## Make a dataframe with the ticks for the chromosome lines
  if (class(region)=="logical") {
    label <- labels$ticks[2]-labels$ticks[1]
    ticks <- numeric()
    chrom <- character()
    y <- numeric()
    for (chr in levels(genome$chrom)) {
      row <- which(genome$chrom == chr)
      tmpTicks <- seq(genome[row,2],genome[row,3],by=label)
      tmpChr <- rep(chr,length(tmpTicks))
      tmpY <- rep(max(bed[bed$chrom==chr,"score"])/12,length(tmpTicks))
      ticks <- append(ticks,tmpTicks)
      chrom <- append(chrom,tmpChr)
      y <- append(y,tmpY)
    }
    chrTicks <- data.frame(chrom=factor(chrom,levels=unique(genome$chrom)),ticks=as.numeric(ticks),y=as.numeric(y))
  }
  ##  plotting
  plot<-ggplot(bed)
  plot<-plot+scale_y_continuous(
    name="Reads per bin",
    limits=c(0,NA)
  )  ##  y scale
  plot<-plot+scale_x_continuous(
    breaks=labels$ticks,
    labels=labels$labels,
    name=xAxisName,
    limits=c(theMin-(theMax-theMin)/100,theMax+(theMax-theMin)/100),expand=c(0.01,0.01)
  )  ##  x scale
  plot<-plot+geom_segment(
    aes(x=chromStart,xend=chromEnd,y=theMedian,yend=theMedian),
    data=genome,
    size=.75,
    color="orchid"
  )  ##  Median line
  if (!plottingRegion) {
    plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=0,yend=0),data=genome,size=.7)  ##  Chromosome length line
    plot <- plot+geom_segment(
      aes(x=ticks,xend=ticks,y=0,yend=y),
      data=chrTicks,size=.7,color="gray25") ## Chromosome line ticks
    plot<-plot+facet_grid(chrom ~ .,scales = "free", space = "free_x")  ##  Facet
    plot<-plot+geom_text(
      aes(x=chromEnd+theMax/200,y=midY,label=chrom),
      data=genome,
      angle=270,
      size=fontSize/pointSize,
      vjust=0
    )  ## Chromosome names
  }
  # if ("cen" %in% colnames(genome)) {
  #   plot<-plot+geom_vline(aes(xintercept=cen),data=genome,color="green3")
  # }
  bedName<-levels(bed$name)[1]  ##  Retrieve sample name
  plot<-plot+ggtitle(paste0(bedName," sequencing coverage"))  ##  Plot title
  plot<-plot+geom_point(aes(x=chromStart+(chromEnd-chromStart)/2,y=score),color="blue",size=0.5)  ##  Raw data
  plot<-plot+gplotTheme
  if (plotting) { plot } else { return(plot) }
}
