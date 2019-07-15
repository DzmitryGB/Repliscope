#' A function to scatterplot 'Trep' column of a Trep dataframe
#' plotTrep function plots values in the 'Trep' column of the supplied dataframe as a function of
#' chromosome coordinates. The genome wide median is plotted as a pink line.
#' @param TrepDF A dataframe containing 'chrom','chromStart', 'chromEnd' and 'Trep' columns (dataframe).
#' @param region Only plot for the provided region in the format 'chrI:1000-3000' (string, optional).
#' @param plotting Should the plot object be sent to the default device? (boolean, defaults to TRUE).
#' @keywords Trep scatterplot genomics bioinformatics
#' @import ggplot2
#' @export
#' @examples
#' plotTrep(TrepDF,region="chrI")

plotTrep <- function(TrepDF,region=FALSE,plotting=TRUE) {
  ## the line below is to appease CRAN checks as it doesn't understand data.frame's variables
  chromStart <- chromEnd <- midY <- Trep <- NULL
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
  theMedian <- as.numeric(summary(TrepDF$Trep)["Median"])
  ymax <- max(TrepDF$Trep,na.rm=T)
  ymin <- min(TrepDF$Trep,na.rm=T)
  ymid <- ymin + (ymax-ymin)/2
  ylims <- c(ymin,ymax+ymax*0.05)
  ##  Create genome dataframe using a helper function
  if (region!=F) {
    plottingRegion <- T
    region <- as.character(region)
    genome <- makeGenome(TrepDF,region)
    if (length(grep(":",region))==1) {
      region <- unlist(strsplit(region,":"))
      tmp <- unlist(strsplit(region[2],"-"))
      region <- c(region[1],tmp)
      TrepDF<-TrepDF[TrepDF$chrom==region[1] & TrepDF$chromStart>=as.integer(region[2]) & TrepDF$chromEnd<=as.integer(region[3]),]
    } else {
      TrepDF<-TrepDF[TrepDF$chrom==region[1],]
    }
  } else {
    genome <- makeGenome(TrepDF)
    plottingRegion <- F
  }
  xAxisName <- attributes(genome)$axisName
  ##  Make nice x-axis labels using a helper function
  theMax <- max(genome$chromEnd)
  theMin <- min(genome$chromStart)
  labels <- makeLabels(theMin,theMax,"b")
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
      tmpY <- rep(max(TrepDF[TrepDF$chrom==chr,"Trep"],na.rm=T)/12,length(tmpTicks))
      ticks <- append(ticks,tmpTicks)
      chrom <- append(chrom,tmpChr)
      y <- append(y,tmpY)
    }
    chrTicks <- data.frame(chrom=factor(chrom,levels=unique(genome$chrom)),ticks=as.numeric(ticks),y=as.numeric(y))
  }
  ##  plotting
  plot<-ggplot(TrepDF)
  plot<-plot+scale_y_reverse(
    name="Trep, min",
    limits=rev(ylims),expand=c(0.1,0.1)
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
    plot<-plot+geom_segment(aes(x=0,xend=chromEnd,y=ylims[2],yend=ylims[2]),data=genome,size=.7)  ##  Chromosome length line
    plot <- plot+geom_segment(
      aes(x=ticks,xend=ticks,y=ylims[2],yend=ymax),
      data=chrTicks,size=.7,color="gray25") ## Chromosome line ticks
    plot<-plot+facet_grid(chrom ~ .,scales = "free_x", space = "free_x")  ##  Facet
    plot<-plot+geom_text(
      aes(x=chromEnd+theMax/200,y=ymid,label=chrom),
      data=genome,
      angle=270,
      size=fontSize/pointSize,
      vjust=0
    )  ## Chromosome names
  }
  plot<-plot+ggtitle("Median replication time (Trep)")  ##  Plot title
  plot<-plot+geom_point(aes(x=chromStart+(chromEnd-chromStart)/2,y=Trep),color="cadetblue",size=0.5)  ##  Raw data
  plot<-plot+gplotTheme
  if (plotting) { plot } else { return(plot) }
}
