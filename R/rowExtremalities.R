#' retrieve row-wise ratios of the actual SDs vs theoretical max (Bernoulli SD)
#' 
#' This is useful in finding probes/features that discriminate classes well.
#'
#' @param x     a Matrix or GenomicRatioSet 
#' 
#' @return      actualSd/bernoulliSd for each row
#'
#' @export
rowExtremalities <- function(x) {
  if (is(x, "GenomicRatioSet")) x <- getBeta(x)
  # set_verbose_block_processing(TRUE) 
  # setAutoBlockSize(1e9) # look at billion entries at a time
  means <- rowMeans2(x, na.rm=TRUE)
  actualSd <- rowSds(x, na.rm=TRUE)
  bernoulliSd <- sqrt(means * (1 - means))
  # practical fix for numerical instability:
  res <- (actualSd / pmax(bernoulliSd, actualSd)) 
  names(res) <- rownames(x)
  return(res)
}
