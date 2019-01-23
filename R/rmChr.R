#' A function to remove single chromosome data from a bed dataframe
#' @param bed A bed dataframe containing 'chrom' column (dataframe, required).
#' @param chr Chromosome to remove (string, required).
#' @keywords outliers BED genomics bioinformatics
#' @export
#' @examples
#' bedDF <- rmChr(W303_S,"chrM") ## removes mitochondria

rmChr <- function(bed,chr) {
  bed<-bed[bed$chrom!=chr,]
  rownames(bed) <- 1:nrow(bed)
  bed$chrom <- factor(bed$chrom, levels = unique(bed$chrom))
  return(bed)
}
