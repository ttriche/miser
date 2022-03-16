#' simple function for plotting a mixture model for sex from XY stats
#' 
#' @param x       a GenomicRatioSet, or an XYstats matrix 
#' @param ...     arguments to pass to Mclust 
#' 
#' @return        an mclust fit object with $sex as fitted sex
#' 
#' @import        mclust
#' 
#' @export
fitSexMix <- function(x, ...) {

  if (is(x, "GenomicRatioSet") | is(x, "SummarizedExperiment")) { 
    if ("XYstats" %in% names(metadata(x))) XY <- as.matrix(metadata(x)$XYstats)
    else XY <- XYstats(x)
    XY <- XY[, colnames(x)]
  } else { 
    XY <- x
  }

  fit <- Mclust(flogit(t(XY[c("XBeta", "YNAfrac"), ])), G=1:2, ...)
  if (fit$G < 2) message("Mixture fit failed; you might only have one sex.")
  fit$sex <- c("M","F")[fit$classification] 
  names(fit$sex) <- rownames(fit$data)
  return(fit)
   
}
