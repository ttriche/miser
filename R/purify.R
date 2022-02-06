#' The world's dumbest deconvolution function
#' 
#' Given x with rows and columns, and purity estimates for each column,
#' `purifyBetas` takes (squeezed) logits of assay(x), divides by `purity`
#' to push values towards the extremes, and then takes (de-squeezed) 
#' expits to put the purified values back on a 0-1 scale. 
#' 
#' `purify` is even dumber and simply divides each assays(x) by purity,
#' unless the argument `betas` is TRUE (the default), in which case it
#' calls `purifyBetas` instead. (This is a methylation-centric package,
#' although VAFs can be treated the same way, hence the default setting.)
#'
#' @param x       a SummarizedExperiment-like object, or a matrix of values
#' @param purity  a vector of estimated purities (`sqz` to 1.00000, inclusive)
#' @param sqz     a squeeze factor to avoid +/-Inf values in `flogit` (1e-6)
#' @param betas   are the values of assay(x) between 0 and 1? (default is TRUE) 
#' 
#' @return        whatever `x` is (matrix or SE), purified (see Details)
#' 
#' @details       If you pass `purity` estimates of 0, you'll get an error.
#'                If `x` is a GenomicRatioSet, only assays(x)$Beta is purified.
#'                No attempt is made to preserve imprinting or any such thing,
#'                although the tendency of imprinted regions to fall near zero
#'                on a logit scale, and the fact that zero times anything is 
#'                still zero, does incidentally encourage this. However, any
#'                resemblance to biological sanity is purely coincidental.
#' 
#' @aliases       purifyBetas 
#' 
#' @import        Matrix
#' @import        minfi
#'
#' @export 
purify <- function(x, purity, sqz=1e-6, betas=TRUE) {

  stopifnot(is(x, "SummarizedExperiment") | is(x, "matrix") | is(x, "Matrix"))
  stopifnot(length(purity) == ncol(x))
  if (betas) { 
    if (is(x, "SummarizedExperiment")) {
      assays(x)$Beta <- .purifyBetas(assays(x)$Beta, purity=purity, sqz=sqz)
    } else { 
      x <- .purifyBetas(x, purity, sqz=sqz)
    }
  } else {
    if (is(x, "SummarizedExperiment")) {
      assays(x) <- lapply(assays(x), purify, purity=purity)
    } else { 
      x <- sweep(x, 2, purity, `/`)
    }
  }
  return(x) 

}


# aliased above 
purifyBetas <- function(x, purity, sqz=1e-6) {
 
  purify(x, purity=purity, sqz=sqz, betas=TRUE) 

}


# helper fn
.purifyBetas <- function(x, purity, sqz=1e-6) {
  
  stopifnot(all(purity >= sqz) & sqz < 1)
  x <- fexpit(.purify(flogit(x, sqz=sqz), purity=purity), sqz=sqz)
  x[x < 0] <- 0
  x[x > 1] <- 1
  return(x) 

}


# helper fn
.purify <- function(x, purity) {
  
  stopifnot(all(purity <= 1) & all(purity > 0))
  sweep(x, 2, purity, `/`)

}
