#' A function to smooth ratio values using cubic smoothing spline
#' smoothRatio function splits values from 'ratio' column by chromosome and based the supplied
#' 'groupMin' and 'split' parameters and then applies smooth.spline() function from R stats package.
#' The supplied dataframe may contain multiple ratios, i.e. ratios produced using multiple replicating
#' samples and/or multiple non-replicating samples. This must be reflected in 'name.rep' and 'name.nonRep'
#' columns. In other words, different ratio dataframes may be merged using rbind() function before calling
#' smoothRatio() function.
#'
#' @param ratioDF A ratio dataframe or combined ratios dataframe containing 'ratio' column (dataframe).
#' @param groupMin Minimum number of values required to make a group (integer, defaults to 5).
#' @param splitNum Minimum number of adjacent bins with missing values to close current group (integer, defaults to 5).
#' @export
#' @importFrom stats na.omit smooth.spline
#' @keywords spline replication genomics bioinformatics
#' @examples
#' ratioDF <- smoothRatio(W303norm)

smoothRatio <- function(ratioDF,groupMin=5,splitNum=5) {
  if (groupMin<4) {
    warning('Minimum group size is four. Setting groupMin parameter to 4.')
    groupMin <- 4
  }
  if (splitNum<1) stop('Split must be a positive integer')
  ratioDF$chrom <- factor(ratioDF$chrom, levels = unique(ratioDF$chrom))
  if ("name.rep" %in% colnames(ratioDF)) {
    ratioDF$name.rep <- factor(ratioDF$name.rep, levels = unique(ratioDF$name.rep))
  }
  if ("name.nonRep" %in% colnames(ratioDF)) {
    ratioDF$name.nonRep <- factor(ratioDF$name.nonRep, levels = unique(ratioDF$name.nonRep))
  }
  ## initiate group number (j)
  j <- 1
  ##  create dummy columns "group" and "splineSmooth"
  if (!("group" %in% colnames(ratioDF))) ratioDF$group <- paste("NA")
  if (!("splineSmooth" %in% colnames(ratioDF))) ratioDF$splineSmooth <- paste("NA")
  ##  initiate vectors
  bin <- vector(mode="numeric",length=0)
  group <- vector(mode="numeric",length=0)
  ##  Actions on every ratio in the supplied dataframe
  for (ratio in unique(paste(as.character(ratioDF$name.rep),as.character(ratioDF$name.nonRep)))) {
    if (length(grep(" ",ratio))==1) {
      ratio <- unlist(strsplit(ratio," "))
      rep <- ratio[1]
      nonRep <- ratio[2]
      #~ 			currentRatio <- ratioDF[ratioDF$name.rep==rep & ratioDF$name.nonRep==nonRep,]
      rows <- which(ratioDF$name.rep==rep & ratioDF$name.nonRep==nonRep)
      chroms <- ratioDF[rows,"chrom"]
      chroms <- factor(chroms,levels=unique(chroms))
      allStarts <- ratioDF[rows,"chromStart"]
      anEnd <- ratioDF[rows,"chromEnd"][1]
      bin <- append(bin,anEnd-allStarts[1])
      if (length(unique(bin))!=1) {
        stop("All samples must have the same genomic bin size")
      } else {
        bin <- unique(bin)
      }
      ##  Actions on every chromosome of the current ratio
      for (chr in levels(chroms)) {
        ## Grouping lines
        starts <- allStarts[chroms==chr]
        group <- j
        i <- 1
        for (i in 1:(length(starts)-1)) {
          i <- i + 1
          if ( (starts[i] - starts[i-1])/bin > (splitNum-1) ) { j <- j+1 }
          group <- append(group,j)
        }
        ## Remove small groups
        for (aGroup in unique(group)) {
          tmp<-which(group == aGroup)
          if (length(tmp) < groupMin) {
            group[tmp] <- NA
          }
        }
        ratioDF[ratioDF$name.rep==rep & ratioDF$name.nonRep==nonRep & ratioDF$chrom==chr,"group"] <- as.integer(group)
        j <- j+1
      }
    }
  }
  #~ Smoothing (spline) in groups
  for (aGroup in as.integer(unique(na.omit(ratioDF$group)))) {
    rows <- which(ratioDF$group==aGroup)
    x <- ratioDF[rows,"chromStart"]
    y <- ratioDF[rows,"ratio"]
    splineObject <- smooth.spline(x,y)
    ratioDF[rows,"splineSmooth"] <- round(as.numeric(splineObject$y),3)
  }
  ratioDF$splineSmooth <- suppressWarnings(as.numeric(ratioDF$splineSmooth))
  ratioDF$group <- suppressWarnings(as.integer(ratioDF$group))
  return(ratioDF)
}
