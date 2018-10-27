#' make sure annotations and so forth are in place before doing a big run...
#' 
#' @param   samps   one or more rows of the `samps` data.frame with `Basename`
#' @param   elts    which elements to extract for the array Basename (1:3)
#'
#' @return          the status of an attempt to annotate those IDATs 
#' 
#' @import  minfi
#' 
#' @export
testIDATs <- function(samps, elts=1:3) {
  rgSet <- read.metharray.exp(base=".", targets=head(samps, 1), verbose=TRUE)
  if (is(sesamize(rgSet), "GenomicRatioSet")) return(TRUE) 
  else stop("Could not map samples to genome.")
} 
