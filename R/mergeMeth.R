#' helper function to merge grSets
#' 
#' @param x   a GenomicRatioSet
#' @param y   another GenomicRatioSet
#' 
#' @return    a GenomicRatioSet combining the two
#' 
#' @import    minfi
#' 
#' @export 
mergeMeth <- function(x, y) { 
  stopifnot(is(x, "GenomicRatioSet"))
  stopifnot(is(y, "GenomicRatioSet"))
  matchingRows <- intersect(rownames(x), rownames(y))
  matchingColData <- intersect(names(colData(x)), names(colData(y)))
  colData(x) <- colData(x)[, matchingColData]
  colData(y) <- colData(y)[, matchingColData]
  cbind(x[matchingRows,], y[matchingRows,])
}
