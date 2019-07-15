#' A helper function to create a gemome dataframe
#'
#' makeGenome is called by plotGenome() and plotCoverage() functions if a genome dataframe is not
#' provided. It creates a BED-like dataframe containing unique chromosome names, their start coordinates
#' (assumed 0), their ends (highest value in the corresponding 'chromEnd' of the BED dataframe) and 'midY'
#' column containing half the max value of the "score" or "ratio" columns per chromosome. This later
#' serves as y coordinate to add chromosome names. Axis name is passed via comment to the output dataframe.
#' Extract it with 'attributes(genome)$axisName'.
#' @param DF A BED or ratio dataframe containing either 'score' or 'ratio' column (dataframe).
#' @param region String in the format 'chrI:1000-3000' (string, optional).
#' @keywords plotting genomics bioinformatics
#' @export
#' @examples
#' genomeDF <- makeGenome(W303_G2)

makeGenome <- function(DF,region=FALSE) {
  ## the line below is to appease CRAN checks as it doesn't understand data.frame's variables
  chrom <- NULL
  if (region==F) {
    ##  using the whole supplied dataframe
    DF$chrom <- factor(DF$chrom,levels=levels(DF$chrom))
    chroms <- levels(DF$chrom)
    ends <- vector(mode="numeric",length=0)
    midY <- vector(mode="numeric",length=0)
    if ("score" %in% colnames(DF)) {
      for (chr in chroms) {
        ends <- append(ends,max(DF[DF$chrom==chr,"chromEnd"]))
        midY <- append(midY,max(DF[DF$chrom==chr,"score"])/2)
      }
    } else if ("ratio" %in% colnames(DF)) {
      for (chr in chroms) {
        ends <- append(ends,max(DF[DF$chrom==chr,"chromEnd"]))
        midY <- append(midY,max(DF[DF$chrom==chr,"ratio"])/2)
      }
    }  else if ("Trep" %in% colnames(DF)) {
      for (chr in chroms) {
        ends <- append(ends,max(DF[DF$chrom==chr,"chromEnd"]))
        midY <- append(midY,max(na.omit(DF[DF$chrom==chr,"Trep"]))/2)
      }
    }
    genome <- data.frame(
      "chrom"=as.character(chroms),
      "chromStart"=as.integer(rep(0,length(chroms))),
      "chromEnd"=as.integer(ends),
      "midY"=as.numeric(midY))
    genome$chrom <- factor(genome$chrom,levels=levels(DF$chrom))
    attr(genome,"axisName") <- paste("Chromosome coordinates")
  } else {
    genome <- data.frame(chrom=as.character(NA),chromStart=as.integer(NA),chromEnd=as.integer(NA),stringsAsFactors=FALSE)
    ##  if region is supplied, first create the region DF
    region <- as.character(region)
    if (length(grep(":",region))==1) {
      region <- unlist(strsplit(region,":"))
      tmp <- unlist(strsplit(region[2],"-"))
      region <- c(region[1],tmp)
    } else {
      tmp <- c(0,max(subset(DF,chrom==region,select="chromEnd")$chromEnd))
      region <- append(region,tmp)
    }
    if (length(region)==3) {
      if (region[3]>region[2]) {
        genome[1,1] <- as.character(region[1])
        genome[1,2] <- as.integer(region[2])
        genome[1,3] <- as.integer(region[3])
      } else {
        stop("Supplied region end should be bigger than region start")
      }
    } else {
      stop("Couldn't understand the supplied region. Please supply it as a string in 'chrX' or 'chrX:1000-2000' format.")
    }
    attr(genome,"axisName") <- as.character(paste(region[1],"coordinates"))
    genome$chrom <- factor(genome$chrom,levels=unique(genome$chrom))
  }
  return(genome)
}
