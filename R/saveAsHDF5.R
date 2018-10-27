#' EXPERIMENTAL harness to eventually bypass in-core processing for HDF5
#' 
#' Right now, all this does is save the GenomicRatioSet to an HDF5 directory.
#' 
#' @param x     a GenomicRatioSet (TODO: support RGChannelSets and SignalSets)
#' @param dir   name of the directory to store the object its data (name of x)
#' @param ...   additional arguments to pass to saveHDF5SummarizedExperiment
#' 
#' @return      invisibly, the same type of object, but now HDF5-backed. 
#'
#' @import      SummarizedExperiment
#' @import      DelayedArray 
#' @import      HDF5Array 
#'
#' @export
saveAsHDF5 <- function(x, dir=NULL, ...) { 
  if (is.null(dir)) dir <- deparse(match.call()$x) 
  saveHDF5SummarizedExperiment(x=x, dir=dir, ...) 
}
