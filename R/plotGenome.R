#' plotGenome: plot replication profile.
#'
#' plotGenome plots scatterplot/barplot/polygon of 'score' and/or 'splineSmooth' columns values
#' by genomic coordinates, either for the whole genome or a region. It also allows annotation
#' with additional BED-formatted data. Ggplot2 package is used for plotting.
#' @param ratioDFs A ratio dataframe or combined ratios dataframe containing 'ratio' column (dataframe).
#' @param geom ggplot2 geom to use for plotting: "geom_point","geom_ribbon" or "geom_segment" (string, defaults to "geom_point").
#' @param ylims A vector of two values for y axis limits - first is lowest, second is highest (numeric vector, defaults to 1 and 2)
#' @param plotting Should the plot object be sent to the default device? (boolean, defaults to TRUE).
#' @param genome A mask dataframe to exclude data from the ratio dataframe (dataframe, optional).
#' The genome dataframe must contain "chrom","chromStart" and "chromEnd" columns.
#' @param region Only plot for the provided region in the format 'chrI:1000-3000' (string, optional).
#' @param guide A dataframe guiding the plotGenome function how to plot the data (dataframe, optional).
#' The guide dataframe must contain the following columns:
#'   'order' (integer) - order to plot data in,
#'   'name.rep' (character) - replicating sample name that matches the one in the ratioDFs dataframe,
#'   'name.nonRep' (character) - non-replicating sample name that matches the one in the ratioDFs dataframe,
#'   'raw' (logical) - should the raw raw data be plotted?
#'   'smooth' (logical) - should the smoothed data be plotted?
#'   'color'(character) - R color to plot the current sample with, both raw and smoothed data.
#' @param lines Additionally plot vertical lines from a BED formatted dataframe (dataframe, optional).
#' @param circles Additionally plot circles on the chromosome line from a BED formatted dataframe (dataframe, optional).
#' @param rectangles Additionally plot rectangles on the chromosome line from a BED formatted dataframe (dataframe, optional).
#' @param pointers Additionally plot downward pointing triangles from a BED formatted dataframe (dataframe, optional).
#' @param colourLines Colour for 'lines' data (string, defaults to green).
#' @param colourCircles Colour for 'circles' data (string, defaults to white).
#' @param colourRectangles Colour for 'rectangles' data (string, defaults to red).
#' @param colourPointers Colour for 'pointers' data (string, defaults to orange).
#' @keywords replication ggplot2 genomics bioinformatics
#' @import ggplot2
#' @importFrom grDevices col2rgb
#' @export
#' @examples
#' plotGenome(sortSeq,geom="geom_ribbon",guide=guide,region="chrIX:250000-439885",
#'     lines=sacCer3[["cen"]],circles=sacCer3[["ori"]])
#' # plot data as polygon for the specified region of chromosome 9
#'
#' plotGenome(syncSeq[["data"]],geom="geom_segment",guide=syncSeq[["guide"]],
#'     region="chrVII:0-1090944",genome=sacCer3[["genome"]],lines=sacCer3[["cen"]],
#'     circles=sacCer3[["ori"]],colourLines="black")
#'
#' plotGenome(MFA,region='chr1:0-2848000')
#' # plot marker frequency analysis for H.volcanii isolate DS2

plotGenome <- function(ratioDFs,geom="geom_point",ylims=c(1,2),plotting=TRUE,
                       genome=NULL,region=FALSE,guide=NULL,
                       lines=NULL,circles=NULL,rectangles=NULL,pointers=NULL,
                       colourLines='#00FF00',colourCircles='#FFFFFF',
                       colourRectangles='#FF0000',colourPointers='#FF7F00'
) {
  ## the line below is to appease CRAN checks as it doesn't understand data.frame's variables
  chromStart <- chromEnd <- NULL

  #~ warning(paste("ratioDFs before:",printColTypes(ratioDFs)))
  ratioDFs$chrom <- factor(ratioDFs$chrom, levels = unique(ratioDFs$chrom))
  ratioDFs$name.rep <- factor(ratioDFs$name.rep, levels=unique(ratioDFs$name.rep))
  ratioDFs$name.nonRep <- factor(ratioDFs$name.nonRep, levels=unique(ratioDFs$name.nonRep))
  #~ warning(paste("ratioDFs after:",printColTypes(ratioDFs)))
  if ("name.rep" %in% colnames(ratioDFs)) {
    ratioDFs$name.rep <- factor(ratioDFs$name.rep, levels = unique(ratioDFs$name.rep))
  }
  if ("name.nonRep" %in% colnames(ratioDFs)) {
    ratioDFs$name.nonRep <- factor(ratioDFs$name.nonRep, levels = unique(ratioDFs$name.nonRep))
  }
  ylims <- ylims[order(ylims)]
  samples <- unique(ratioDFs[,c("name.rep","name.nonRep")])
  rownames(samples) <- 1:nrow(samples)
  samples$name <- paste0(as.character(samples[,1])," (",as.character(samples[,2]),")")
  ##  First ratio
  firstRatioName <- samples[1,3]

  #~ 	if (length(grep(" ",firstRatioName))==1) {
  #~ 		ratio <- unlist(strsplit(firstRatioName," "))
  rep <- samples[1,1]
  nonRep <- samples[1,2]
  firstRatio <- ratioDFs[ratioDFs$name.rep==rep & ratioDFs$name.nonRep==nonRep,]
  firstRatio$chrom <- factor(firstRatio$chrom,levels=levels(ratioDFs$chrom))
  #~  		warning(paste0("First ratio: ",printColTypes(firstRatio)))
  #~ 		warning(paste0("First ratio: ",dfHead(firstRatio)))
  #~ 	}

  ##  If region is supplied, create the wee genome dataframe
  if (class(region)=="character") {
    genome <- makeGenome(ratioDFs,region=region)
    xAxisName <- attributes(genome)$axisName
    # region <- as.character(region)
    # if (length(grep(":",region))==1) {
    region <- unlist(strsplit(unlist(strsplit(region,"-")),":"))
    #   genome <- data.frame(chrom=as.character(region[1]),chromStart=as.integer(region[2]),chromEnd=as.integer(region[3]))
    #   xAxisName <- paste(region[1],"coordinates")
    # }
  } else {
    if (is.null(genome)) {
      ##  Genome is not supplied
      ## create genome dataframe using a helper function based on the first sample
      genome <- makeGenome(ratioDFs)
    } else { ## genome dataframe is supplied
      ##  remove missing chromosomes from the genome
      goods<-intersect(levels(firstRatio$chrom),levels(genome$chrom))
      if (!all(levels(genome$chrom) %in% goods)) {
        genome<-genome[genome$chrom %in% goods,]
        genome$chrom<-factor(genome$chrom,levels=unique(genome$chrom))
        rownames(genome) <- 1:nrow(genome)
      }
    }
    xAxisName <- "Chromosome coordinates"
  }
  dummy <- genome ## dummy dataframe
  dummy$ratio <- NA
  #~ 	warning(paste0("Genome: ",printColTypes(genome)))
  #~ 	warning(paste0("Genome: ",dfHead(genome)))

  ##  Make nice x-axis labels using a helper function
  theMax <- max(genome$chromEnd)
  theMin <- min(genome$chromStart)
  labels <- makeLabels(theMin,theMax,"b")

  ## Make a dataframe with the ticks for the chromosome lines
  if (class(region)=="logical") {
    label <- labels$ticks[2]-labels$ticks[1]
    ticks <- numeric()
    chrom <- character()
    for (chr in levels(genome$chrom)) {
      row <- which(genome$chrom == chr)
      tmpTicks <- seq(genome[row,2],genome[row,3],by=label)
      tmpChr <- rep(chr,length(tmpTicks))
      ticks <- append(ticks,tmpTicks)
      chrom <- append(chrom,tmpChr)
    }
    chrTicks <- data.frame(chrom=factor(chrom,levels=unique(genome$chrom)),ticks=as.numeric(ticks))
  }
  ##  Get the bin size for the first sample
  bin <- firstRatio$chromEnd[1]-firstRatio$chromStart[1]

  ##  load ggplot2 library and set plotting theme
  fontSize <- 10
  pointSize <- 1 / 0.352777778
  gplotTheme<-theme(
    text = element_text(size=fontSize),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5,size=rel(1.8)),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour="gray50",size=.5),
    axis.title = element_text(size = rel(1.6),face="bold"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=rel(1.4)),
    legend.position="bottom",
    legend.spacing=unit(0.1,"cm"),
    legend.title=element_blank(),
    legend.text = element_text(size=rel(1)),
    legend.key = element_blank(),
    legend.key.size = unit(c(0.4,0.4),"cm"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")
  )

  ##  start creating plot object - labels, scales names etc.
  plot <- ggplot(genome)
  plot <- plot + scale_y_continuous(
    name="Relative copy number",
    limits = c(ylims[1]-((ylims[2]-ylims[1])/5),ylims[2]+((ylims[2]-ylims[1])/5))
  )  ##  y scale
  plot <- plot + annotate(
    geom="segment",
    x=theMin,
    xend=theMin,
    y=ylims[1],
    yend=ylims[2],
    colour="gray30",
    lwd=0.5
  )  ##  y scale line
  plot <- plot + scale_x_continuous(
    breaks=labels$ticks,
    labels=labels$labels,
    name=xAxisName,
    limits=c(theMin-(theMax-theMin)/50,theMax+(theMax-theMin)/50),
    expand=c(0.02,0.02)
  )  ##  x scale
  plot <- plot + geom_segment(
    aes(x=theMin,xend=genome$chromEnd,y=(ylims[1]-((ylims[2]-ylims[1])/10)),yend=(ylims[1]-((ylims[2]-ylims[1])/10))),
    size=.8,color="gray20"
  )  ##  Chromosome length line
  if (class(region)=="logical") {
    plot <- plot + facet_grid(
      chrom ~ .,
      scales = "free",
      space = "free_x"
    )  ##  Facet
    plot <- plot+geom_text(
      aes(x=chromEnd+(theMax-theMin)/100,y=ylims[1]+((ylims[2]-ylims[1])/2),label=chrom),
      angle=270,
      size=fontSize/pointSize,
      vjust=0
    )  ##  Chromosome labels
    plot <- plot+geom_segment(
      aes(x=ticks,xend=ticks,y=ylims[1]-((ylims[2]-ylims[1])/5),yend=ylims[1]-((ylims[2]-ylims[1])/10)),
      data=chrTicks,size=.7,color="gray20") ## Chromosome line ticks
  }
  plot <- plot + geom_segment(
    aes(x=theMin,xend=chromEnd,y=ylims[1],yend=ylims[1]),
    linetype="dashed",
    size=.25,
    color="gray20"
  )  ##  Lower y line
  plot <- plot + geom_segment(
    aes(x=theMin,xend=chromEnd,y=ylims[2],yend=ylims[2]),
    linetype="dashed",
    size=.25,
    color="gray20"
  )  ##  Upper y line
  plot <- plot + geom_text(
    aes(x=theMin-(theMax-theMin)/100,y=ylims[1],label=as.character(ylims[1])),
    size=fontSize/pointSize,
    hjust=1
  )  ##  Lower y break
  plot <- plot + geom_text(
    aes(x=theMin-(theMax-theMin)/100,y=ylims[2],label=as.character(ylims[2])),
    size=fontSize/pointSize,
    hjust=1
  )  ##  Upper y break

  ##  initialising legendary vectors
  color.names <- vector(mode="character",length=0)
  color.values <- vector(mode="character",length=0)
  color.lines <- vector(mode="character",length=0)
  color.shapes <- vector(mode="numeric",length=0)
  color.sizes <- vector(mode="numeric",length=0)
  color.fills<-  vector(mode="character",length=0)


##  unguided plotting: create guide DF
  if (is.null(guide)) {
	##  Colors
    dotColors <- c("gray50","deepskyblue4","orange3","darkolivegreen4","hotpink4","purple4","red3","yellow3")
    ##  Check for smoothed data
    if ("splineSmooth" %in% colnames(ratioDFs)) {
  		smooth <- TRUE
      } else {
  		smooth <- FALSE
  	}
    ##  Guide DF
  	guide <- data.frame(
  		order=as.integer(1:nrow(samples)),
  		name.rep=as.character(samples$name.rep),
  		name.nonRep=as.character(samples$name.nonRep),
  		raw=as.logical(rep(TRUE,nrow(samples))),
  		smooth=as.logical(rep(smooth,nrow(samples))),
  		color=as.character(dotColors[1:nrow(samples)]),
  		stringsAsFactors=F)
    }

## guided plotting
	geoms <- c("geom_point","geom_ribbon","geom_segment")
	geom_pars <- c(
	  "geom_point(aes(x=chromStart+bin/2,y=ratio,color='changeMe'),size=0.2,data=currentRatio)",
	  "geom_ribbon(aes(x=chromStart+bin/2,ymin=ylims[1],ymax=ratio,color='changeMe'),fill=guide$color[i],data=currentRatio)",
	  "geom_segment(aes(x=chromStart+bin/2,xend=chromStart+bin/2,yend=ratio,y=ylims[1]),data=currentRatio,color=guide$color[i])"
	)
	geom_str <- geom_pars[match(geom,geoms)]
    guide <- guide[order(guide$order,decreasing=T,na.last = T),]
    rownames(guide) <- 1:nrow(guide)
    layers <- length(na.omit(guide$order))
    sampleNames <- character()
    for (i in 1:layers) {
      if (guide[i,"raw"] | guide[i,"smooth"]) {
        rep <- guide$name.rep[i]
        nonRep <- guide$name.nonRep[i]
        currentRatio <- ratioDFs[ratioDFs$name.rep==rep & ratioDFs$name.nonRep==nonRep,]
        ratioFactor <- currentRatio$ratioFactor[1]
        sampleName <- paste0(rep," (",nonRep,",",ratioFactor,")")
        sampleNames <- append(sampleNames,sampleName)
        if (class(region)=="character") {
          currentRatio <- currentRatio[currentRatio$chrom==region[1] & currentRatio$chromStart>=as.numeric(region[2]) & currentRatio$chromEnd<=as.numeric(region[3]),]
        }
        if (xor(as.logical(guide$raw[i]),as.logical(guide$smooth[i]))) { ##  Either raw only or smooth only
          geom_string <- gsub("changeMe",sampleName,geom_str)
          if (as.logical(guide$smooth[i])) { # smooth only
            geom_string <- gsub("=ratio,","=splineSmooth,group=group,",geom_string)
            geom_string <- gsub("geom_point","geom_line",geom_string)
            geom_string <- gsub("size=0.2,","",geom_string)
          }
          plot <- plot + eval(parse(text=geom_string))
        } else if (as.logical(guide$raw[i]) & as.logical(guide$smooth[i])) { ## Both
          geom_string <- gsub("changeMe",sampleName,geom_str)
          plot <- plot + eval(parse(text=geom_string))

          ##  make smoothed line darker
          rgbColor <- col2rgb(guide$color[i])
          lineColor <- rgbColor-25
          hexVector <- character()
          for (j in 1:3) {
            if (is.na(suppressWarnings(sqrt(lineColor[j])))) {
              lineColor[j] <- 0
              hexVector <- append(hexVector,"00")
            } else {
              hexByte <- as.character(as.hexmode(lineColor[j]))
              if (nchar(hexByte)==1) {
                hexByte <- paste0("0",hexByte)
              }
              hexVector <- append(hexVector,hexByte)
            }
          }
          lineColorHex <- paste0("#",paste0(hexVector,collapse=""))
          ##  add smoothed line to the plot object
          geom_string <- paste0("geom_line(aes(x=chromStart+bin/2,y=splineSmooth,group=group),data=currentRatio,color='",lineColorHex,"')")
          plot <- plot + eval(parse(text=geom_string))
        }
        color.names<-c(color.names,sampleName)
        if (geom == "geom_point") {
          color.values<-c(color.values,as.character(guide$color[i]))
          color.shapes<-c(color.shapes,20)
          color.lines<-c(color.lines,"solid")
          color.sizes<-c(color.sizes,5)
          color.fills<-c(color.fills,NA)
        } else if (geom == "geom_ribbon") {
			color.fills<-c(color.fills,as.character(guide$color[i]))
			color.shapes<-c(color.shapes,22)
			color.sizes<-c(color.sizes,4)
			color.values<-c(color.values,NA)
			color.lines<-c(color.lines,"blank")
        } else if (geom=="geom_segment") {
          geom_string <- paste0("geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,color='",sampleName,"'),data=dummy)")
          plot <- plot + eval(parse(text=geom_string))
			color.lines<-c(color.lines,"solid")
			color.values<-c(color.values,as.character(guide$color[i]))
			color.shapes<-c(color.shapes,NA)
			color.sizes<-c(color.sizes,2)
			color.fills<-c(color.fills,NA)
        }
      }
    }

  ##  optional plot features
  if ("p.value" %in% colnames(firstRatio)) {
    temp<-which(firstRatio$p.value < 0.001)
    sign.999<-firstRatio[temp,]
   # sign.999$chrom <- factor(sign.999$chrom,levels=levels(ratioDFs))
    temp<-which(firstRatio$p.value <= 0.01 & firstRatio$p.value > 0.001)
    sign.99<-firstRatio[temp,]
    ##  Plot significant p-values
    if (nrow(sign.999)!=0) {
      if (class(region)=="character") {
        sign.999 <- sign.999[sign.999$chrom==region[1] & sign.999$chromStart>=as.numeric(region[2]) & sign.999$chromEnd<=as.numeric(region[3]),]
      }
      plot<-plot+geom_segment(aes(
        x=chromStart+bin/2,
        xend=chromStart+bin/2,
        y=ylims[2]+((ylims[2]-ylims[1])/10),
        yend=ylims[2]+((ylims[2]-ylims[1])/5)
      ),data=sign.999,color="gray25")
      plot<-plot+geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,color="p-value***"),data=dummy,show.legend=TRUE)  ##  dummy
      #~ 			colors <- rbind(colors,data.frame(name="p-value***",color="gray25",shape=NA,line="solid",size=1,stringsAsFactors=F))
      color.names<-c(color.names,"p-value***")
      color.lines<-c(color.lines,"solid")
      color.values<-c(color.values,"gray25")
      color.shapes<-c(color.shapes,NA)
      color.sizes<-c(color.sizes,2)
      color.fills<-c(color.fills,NA)
    }
    if (nrow(sign.99)!=0) {
      if (class(region)=="character") {
        sign.99 <- sign.99[sign.99$chrom==region[1] & sign.99$chromStart>=as.numeric(region[2]) & sign.99$chromEnd<=as.numeric(region[3]),]
      }
      plot<-plot+geom_segment(aes(
        x=chromStart+bin/2,
        xend=chromStart+bin/2,
        y=ylims[2]+((ylims[2]-ylims[1])/10),
        yend=ylims[2]+((ylims[2]-ylims[1])/5)
      ),data=sign.99,color="gray75")
      plot<-plot+geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,color="p-value**"),data=dummy,show.legend=TRUE)  ##  dummy
      #~ 			colors <- rbind(colors,data.frame(name="p-value**",color="gray50",shape=NA,line="solid",size=1,stringsAsFactors=F))
      color.names<-c(color.names,"p-value**")
      color.lines<-c(color.lines,"solid")
      color.values<-c(color.values,"gray75")
      color.shapes<-c(color.shapes,NA)
      color.sizes<-c(color.sizes,2)
      color.fills<-c(color.fills,NA)
    }
  }

  if (!is.null(lines)) {
    lines$chrom <- factor(lines$chrom,levels=levels(genome$chrom))
    linesName <- as.character(lines$name[1])
    if (class(region)=="character") {
      lines <- lines[lines$chrom==region[1] & lines$chromStart>=as.numeric(region[2]) & lines$chromEnd<=as.numeric(region[3]),]
    }
    plot<-plot+geom_segment(aes(x=chromStart,xend=chromEnd,y=ylims[1],yend=ylims[2]),data=lines,color=colourLines)
    plot<-plot+geom_vline(aes(xintercept=theMin-(theMax-theMin)/10,color=linesName),data=dummy,show.legend=TRUE)  ##  dummy
    color.names<-c(color.names,linesName)
    color.lines<-c(color.lines,"solid")
    color.values<-c(color.values,colourLines)
    color.shapes<-c(color.shapes,NA)
    color.sizes<-c(color.sizes,2)
    color.fills<-c(color.fills,NA)
  }

  if (!is.null(circles)) {
    circles$chrom <- factor(circles$chrom,levels=levels(genome$chrom))
    circlesName <- as.character(circles$name[1])
    if (class(region)=="character") {
      circles <- circles[circles$chrom==region[1] & circles$chromStart>=as.numeric(region[2]) & circles$chromEnd<=as.numeric(region[3]),]
    }
    plot <- plot+geom_point(aes(x=chromStart+(chromEnd-chromStart)/2,y=ylims[1]-((ylims[2]-ylims[1])/10),color=circlesName),
                            data=circles,fill=colourCircles,size=2,shape=21)
    color.names<-c(color.names,circlesName)  ##  Circles legend values
    color.fills<-c(color.fills,colourCircles)
    color.shapes<-c(color.shapes,21)
    color.sizes<-c(color.sizes,3)
    color.values<-c(color.values,"black")
    color.lines<-c(color.lines,"blank")
  }

  if (!is.null(rectangles)) {  ##  If a dataframe supplied, plot the rectangles
    rectangles$chrom <- factor(rectangles$chrom,levels=levels(genome$chrom))
    rectanglesName <- levels(rectangles$name)[1]
    if (class(region)=="character") {
      rectangles <- rectangles[rectangles$chrom==region[1] & rectangles$chromStart>=as.numeric(region[2]) & rectangles$chromEnd<=as.numeric(region[3]),]
    }
    plot <- plot+geom_rect(aes(
      xmin=chromStart,xmax=chromEnd,
      ymin=(ylims[1]-((ylims[2]-ylims[1])/10)) - (ylims[2]-ylims[1])/20,
      ymax=(ylims[1]-((ylims[2]-ylims[1])/10)) + (ylims[2]-ylims[1])/20),
      data=rectangles,color=colourRectangles,fill=NA)
    plot <- plot+geom_point(aes(x=theMin-(theMax-theMin)/10,y=ylims[1]-((ylims[2]-ylims[1])/2),color=rectanglesName),data=dummy,fill=NA,shape=22)  ##  dummy
    color.names<-c(color.names,rectanglesName)  ##  Rectangles legend values
    color.fills<-c(color.fills,NA)
    color.shapes<-c(color.shapes,22)
    color.sizes<-c(color.sizes,5)
    color.values<-c(color.values,colourRectangles)
    color.lines<-c(color.lines,"blank")
  }
  if (!is.null(pointers)) {  ##  If a dataframe supplied, plot the pointers
    pointers$chrom <- factor(pointers$chrom,levels=levels(genome$chrom))
    pointersName <- levels(pointers$name)[1]
    if (class(region)=="character") {
      pointers <- pointers[pointers$chrom==region[1] & pointers$chromStart>=as.numeric(region[2]) & pointers$chromEnd<=as.numeric(region[3]),]
    }
    plot <- plot+geom_point(aes(x=chromStart+(chromEnd-chromStart)/2,y=ylims[1],color=pointersName),
                            data=pointers,fill=colourPointers,size=3,shape=25)
    color.names<-c(color.names,pointersName)  ##  Pointers legend values
    color.fills<-c(color.fills,colourPointers)
    color.shapes<-c(color.shapes,25)
    color.sizes<-c(color.sizes,4)
    color.values<-c(color.values,"black")
    color.lines<-c(color.lines,"blank")
  }

  my.colors<-color.values  ##  Preparing named vector for colors
  names(my.colors)<-color.names

#~   plot <- plot + scale_colour_manual(
#~     name="",values=rev(my.colors),guide = guide_legend(override.aes = list(
#~ 	  linetype=rev(color.lines),
#~ 	  shape=rev(color.shapes),
#~ 	  size=rev(color.sizes),
#~ 	  fill=rev(color.fills)),reverse=T,nrow=2),breaks=names(my.colors)
#~   )

  plot <- plot + scale_colour_manual(name='legend',values=rev(my.colors),guide=
		guide_legend(override.aes = list(
#~ 						values=rev(my.colors),
						linetype=rev(color.lines),
						shape=rev(color.shapes),
						size=rev(color.sizes),
						fill=rev(color.fills)
		),reverse=T,nrow=2),breaks=names(my.colors))
  plot <- plot + scale_fill_manual(name='legend',values=rev(color.fills),guide='none')
  plot <- plot + scale_shape_manual(name='legend',values=rev(color.shapes),guide='none')
  plot <- plot + scale_size_manual(name='legend',values=rev(color.sizes),guide='none')
  plot <- plot + scale_linetype_manual(name='legend',values=rev(color.lines),guide='none')

#~   plot <- plot + guides(
#~ 		shape=guide_legend(override.aes = list(
#~ 			values=rev(color.shapes)
#~ 		)),
#~ 		fill=guide_legend(override.aes = list(
#~ 			values=rev(color.fills)
#~ 		)),
#~ 		linetype=guide_legend(override.aes = list(
#~ 			values=rev(color.lines)
#~ 		)),
#~ 		size=guide_legend(override.aes = list(
#~ 			values=rev(color.sizes)
#~ 		)),
#~ 		colour = guide_legend(override.aes = list(
#~ 				values=rev(my.colors)
#~ 				linetype=rev(color.lines),
#~ 				shape=rev(color.shapes),
#~ 				size=rev(color.sizes),
#~ 				fill=rev(color.fills)
#~ 			))
#~ 		)

  ##  custom theme
  plot <- plot + gplotTheme

  ##  and we're done
  if (plotting) {
#~    if (region==F) {
#~      X11(width=7,height=10,pointsize=15)
#~    } else {
#~      X11(width=7,height=3,pointsize=15)
#~    }
    plot
  } else {
    return(plot)
  }
}
