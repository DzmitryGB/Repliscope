#' A function to calculate 'score' ratio between two bed dataframes
#' makeRatio merges two supplied bed dataframes, calculates ratio of their "score" values
#' normalises the ratio by the 'score' sums.
#' @param bedRep Bed dataframe containing read counts from a replicating sample (dataframe).
#' The bed dataframe must contain "chrom","chromStart", "chromEnd" and "score" columns.
#' @param bedNonRep Bed dataframe containing read counts from a non-replicating sample (dataframe).
#' The bed dataframe must contain "chrom","chromStart", "chromEnd" and "score" columns.
#' @keywords replication genomics bioinformatics
#' @export
#' @examples
#' ratioDF <- makeRatio(W303_S,W303_G2)

makeRatio <- function(bedRep,bedNonRep) {
  # get intersected chromosomes for later
  chromRep <- unique(bedRep$chrom)
  chromNonRep <- unique(bedNonRep$chrom)
  chroms <- intersect(chromRep,chromNonRep)
  byVector <- c("chrom","chromStart","chromEnd")
  ratioDF <- merge(bedRep,bedNonRep,by=byVector,all=F,sort=F)
  # rename columns
  mergedNames <- (colnames(ratioDF))
  mergedNames <- gsub(".x",".rep",mergedNames)
  mergedNames <- gsub(".y",".nonRep",mergedNames)
  names(ratioDF) <- mergedNames
  #~ Calculate ratio
  repSum <- sum(as.numeric(ratioDF$score.rep))
  nonRepSum <- sum(as.numeric(ratioDF$score.nonRep))
  corrFactor <- repSum/nonRepSum
  ratioDF$ratio <- round(ratioDF$score.rep/ratioDF$score.nonRep/corrFactor,4)
  #~ Clean up the dataframe
  ratioDF <- ratioDF[,!(names(ratioDF) %in% c("score.rep","score.nonRep"))]
  ## order the resulting ratio dataframe
  # ratioDF<-ratioDF[order(ratioDF$chrom,ratioDF$chromStart),]
  # rownames(ratioDF) <- 1:nrow(ratioDF)
  ratioDF$ratioFactor <- "1.0"
  ratioDF$ratioFactor <- factor(ratioDF$ratioFactor,levels="1.0")
  ratioDF$chrom <- factor(ratioDF$chrom,levels=chroms)
  return(ratioDF)
}
