#' helper function to merge grSets while keeping their SNP matrices intact
#' 
#' @param x   a GenomicRatioSet
#' @param y   another GenomicRatioSet
#' 
#' @return    a GenomicRatioSet with sensibly managed SNPs in its metadata()
#' 
#' @import minfi
#' 
#' @export 
mergeWithSNPs <- function(x, y) { 
  stopifnot(is(x, "GenomicRatioSet"))
  stopifnot(is(y, "GenomicRatioSet"))
  z <- mergeMeth(x, y)
  metadata(z) <- list() 
  xSNP <- metadata(x)$SNPs[, colnames(x)] 
  ySNP <- metadata(y)$SNPs[, colnames(y)] 
  matchingSNPs <- intersect(rownames(xSNP), rownames(ySNP))
  metadata(z)$SNPs <- cbind(xSNP[matchingSNPs,], ySNP[matchingSNPs,])
  return(z) 
}
