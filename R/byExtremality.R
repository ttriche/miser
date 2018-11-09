#' For array-based DNA methylation, particularly summarized across regions,
#' we can do better (a lot better) than MAD.  Since we know that 
#'
#' max(SD(X_j)) if X_j ~ Beta(a, b) < max(SD(X_j)) if X_j ~ Bernoulli(a/(a+b))
#'
#' for X having a known mean and SD, hence solvable for a + b by MoM, define
#'
#' extremality = sd(X_j) / bernoulliSD(mean(X_j))
#'
#' This function selects the k most extremal rows of x and returns their values.
#'
#' TODO: operate in parallel when working on an HDF5-backed object
#'
#' @param     x   a matrix of beta values, or a GenomicRatioSet
#' @param     k   how many rows to return (500)
#' 
#' @return    the most extremal _k_ rows of _x_ (returns the same class as _x_)
#' 
#' @import    matrixStats
#' @import    DelayedMatrixStats
#' 
#' @export
byExtremality <- function(x, k=500) {
  k <- min(nrow(x), k)
  extremality <- .extremality(x)
  x[rev(order(extremality))[seq_len(k)], ]
}

# helper fn
.extremality <- function(x) {
  if (is(x, "GenomicRatioSet")) {
    DelayedArray:::set_verbose_block_processing(TRUE) 
    setAutoBlockSize(1e6) # look at million entries at a time
    .extremality(getBeta(x))
  } else { 
    means <- DelayedMatrixStats::rowMeans2(x, na.rm=TRUE)
    bernoulliSd <- sqrt(means * (1 - means))
    actualSd <- DelayedMatrixStats::rowSds(x, na.rm=TRUE)
    return(actualSd / bernoulliSd)
  }
}
