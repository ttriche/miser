#' make sure annotations and so forth are in place before doing a big run...
#' 
#' @param   samps   two or more rows of the `samps` data.frame with `Basename`
#' @param   elts    which elements to extract for the array Basename (1:3)
#'
#' @return          the status of an attempt to annotate those IDATs 
#' 
#' @import  minfi
#' 
#' @export
testIDATs <- function(samps, elts=1:3) {
  rgSet <- read.metharray.exp(base=".", targets=samps[1:2,], verbose=TRUE)
  colnames(rgSet) <- rgSet$subject
  res <- sesamize(rgSet)
  is(res, "GenomicRatioSet")
} 
