#' EXPERIMENTAL harness to eventually bypass in-core processing for HDF5
#' 
#' Right now, all this does is load the GenomicRatioSet in an HDF5 directory.
#' 
#' @param dir   name of the directory where the object and its data are stored
#' 
#' @return      HDF5-backed GenomicRatioSet (or SummarizedExperiment, or w/e)
#'
#' @import      SummarizedExperiment
#' @import      DelayedArray 
#' @import      HDF5Array 
#'
#' @export
loadAsHDF5 <- function(dir) loadHDF5SummarizedExperiment(dir=dir)
