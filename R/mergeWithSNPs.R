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
  xSNP <- metadata(x)$SNPs[, colnames(x)] 
  ySNP <- metadata(y)$SNPs[, colnames(y)] 
  matchingRows <- intersect(rownames(x), rownames(y))
  matchingSNPs <- intersect(rownames(xSNP), rownames(ySNP))
  matchingColData <- intersect(names(colData(x)), names(colData(y)))
  colData(x) <- colData(x)[, matchingColData]
  colData(y) <- colData(y)[, matchingColData]
  z <- cbind(x[matchingRows,], y[matchingRows,])
  metadata(z) <- NULL 
  metadata(z)$SNPs <- cbind(xSNP[matchingSNPs,], ySNP[matchingSNPs,])
  return(z) 
}
