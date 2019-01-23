#' Load a BED formatted file.
#'
#' The function reads BED formatted files. The BED file format defined by UCSC: http://genome.ucsc.edu/FAQ/FAQformat.
#'  First three columns ("chrom", "chromStart", "chromEnd") are mandatory. The file fields may be separated by tabs,
#'   spaces or commas. If the BED file contains a header, it will be ignored. If a genome mask dataframe is provided,
#'   only data intersected with the mask will be retained. Resulting data is ordered by "chromStart" columns.
#' @param file Path to the BED file (string, mandatory)
#' @param genome A mask dataframe to exclude data from the BED file (dataframe, optional).
#' The genome dataframe must contain "chrom" column and may further contain "chromStart" and "chromEnd" columns in this
#'  order.
#' @param name A string to replace the 'name' column of the loaded BED file with (string, optional).
#' @keywords BED genomics bioinformatics
#' @export
#' @importFrom utils read.table
#' @examples
#' W303_G2 <- loadBed(system.file("extdata/W303_G2.bed",package="Repliscope"), name='W303_G2')
#' W303_G2_chrI <- loadBed(system.file("extdata/W303_G2.bed",package="Repliscope"),
#'                        name='W303_G2',genome=sacCer3[["genome"]])

loadBed <- function(file,genome=NULL,name=NULL) {
  ## the line below is to appease CRAN checks as it doesn't understand data.frame's variables
  chrom <- NULL
  colNames <- c("chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
  colTypes <- c("factor","integer","integer","factor","integer","factor","integer","integer","string","integer","string","string")
  fileConn <- file(file,'r')
  firstLine <- readLines(fileConn,n=1)
  close(fileConn)
  firstLineVector <- unlist(strsplit(firstLine,c("\t"," ")))
  colNumber <- length(firstLineVector)
  if (colNumber<3) stop("BED file must contain at least 3 columns. Valid field separators are tab and space.")
  if (is.na(suppressWarnings(as.numeric(firstLineVector[2])))) {
    bed<-read.table(file,sep="",skip=1,stringsAsFactors=F,colClasses=colTypes[1:colNumber])
  } else {
    bed<-read.table(file,sep="",stringsAsFactors=F,colClasses=colTypes[1:colNumber])
  }
  colnames(bed)<-colNames[1:colNumber]
  ## replace content of the 'name' column if $name is supplied
  if (!is.null(name)) {
    bed$name<-paste(name)
    bed$name<-factor(bed$name,levels=name)
  }
  ##  original chromosome order
  chroms <- unique(bed$chrom)
  ## if a mask genome dataframe is provided, remove data omitted in it
  if (!is.null(genome)) {
    chroms <- as.character(genome[,1])
    ##  remove whole chromosomes
    junk<-setdiff(levels(bed$chrom),chroms)
    for (bad.chr in junk) {
      trash<-which(bed$chrom == bad.chr)
      if (length(trash) > 0) {
        bed<-bed[-trash,]
      }
    }
    if (length(chroms)>=length(junk)) {
      message("\n\tExcluded data for ",paste0(junk,collapse=","),".\n")
    } else {
      message("\n\tLoaded data for ",paste0(chroms,collapse=","),".\n")
    }
        ##  remove data based on supplied chromosome starts and ends
    if(ncol(genome)>1) {
      for (chr in chroms) {
        chromStart <- as.integer(genome[genome[,1]==chr,2])
        chromEnd <- as.integer(genome[genome[,1]==chr,3])
        bed <- bed[!(bed$chrom==chr & (bed$chromStart<chromStart | bed$chromEnd>chromEnd)),]
      }
    }
  }
  ## order the resulting bed dataframe
  orderedBed <- bed[0,]
  for (chr in chroms) {
    tmp <- subset(bed,chrom==chr)
    orderedBed <- rbind(orderedBed,tmp[order(tmp$chromStart),])
  }
  orderedBed$chrom <- factor(orderedBed$chrom, levels=chroms)
  rownames(orderedBed) <- 1:nrow(orderedBed)
  return(orderedBed)
}
