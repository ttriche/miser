#' retrieve row-wise ratios of the actual SDs vs theoretical max (Bernoulli SD)
#' 
#' This is useful in finding probes/features that discriminate classes well.
#'
#' @param x     a Matrix, DelayedMatrix, or GenomicRatioSet 
#' 
#' @return      actualSd/bernoulliSd for each row
#'
#' @import DelayedMatrixStats
#' @import DelayedArray
#'
#' @export
rowExtremalities <- function(x) {
  if (is(x, "GenomicRatioSet")) {
    x <- getBeta(x)
  } 
  DelayedArray:::set_verbose_block_processing(TRUE) 
  setAutoBlockSize(1e6) # look at million entries at a time
  means <- DelayedMatrixStats::rowMeans2(x, na.rm=TRUE)
  actualSd <- DelayedMatrixStats::rowSds(x, na.rm=TRUE)
  bernoulliSd <- sqrt(means * (1 - means))
  # practical fix for numerical instability:
  return(actualSd / pmax(bernoulliSd, actualSd)) 
}
