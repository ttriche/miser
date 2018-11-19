#' The world's dumbest deconvolution function
#' 
#' Given estimated purities for each sample, and a GenomicRatioSet,
#' `purify` takes (squeezed) logits of beta values, divides by purity 
#' to push values towards the extremes, and then takes (de-squeezed) 
#' expits to put the purified values back on a 0-1 scale. No attempt 
#' is made to preserve imprinting or any such thing, although the 
#' tendency of imprinted regions to end up near zero on a logit scale
#' and the tendency of zeroes multiplied by anything to be zeroes does
#' incidentally encourage this. However, any resemblance to biological 
#' sanity is purely coincidental. Also, we don't care how you estimate 
#' purity. Spin the bottle, use flow cytometry, whatever you like. All
#' this function does is take logits, divide, and take expits.
#' 
#' @param x       a GenomicRatioSet (if passed a matrix, will treat as betas)
#' @param purity  a vector of estimated purities (`sqz` to 1.00000, inclusive)
#' @param sqz     a squeeze factor to avoid +/-Inf values in `flogit` (1e-6)
#' 
#' @return        whatever `x` is (matrix or grSet), with purified beta values 
#' 
#' @import        minfi
#'
#' @export 
purifyBetas <- function(x, purity,  sqz=0.000001) { 

  stopifnot(all(purity >= sqz) & all(purity <= 1))
  stopifnot(length(purity) == ncol(x))

  if (is(x, "GenomicRatioSet")) {
    message("Updating beta values to purified estimates...") 
    assays(x)$Beta <- .purify(getBeta(x), purity, sqz=sqz)
  } else if (is(x, "matrix")) { 
    message("Estimating purified beta values...") 
    x <- .purify(x, purity, sqz=sqz)
  }
  message("...done.") 
  return(x) 

}


# helper fn
.purify <- function(x, purity, sqz=0.000001) {
  x <- t(fexpit(t(flogit(x, sqz=sqz))/purity, sqz=sqz)) 
  x[x < 0] <- 0
  x[x > 1] <- 1 
  return(x) 
}
