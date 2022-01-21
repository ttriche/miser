#' read in a directory's worth of IDATs, flag controls, and sesamize them all
#' 
#' Note: this function simply reads all IDATs in the current directory, 
#'       runs QC, stuffs that into metadata, and sesamizes the result. 
#'
#' @param ...   options to pass to sesame::sesamize
#'
#' @return a GenomicRatioSet 
#' 
#' @export 
processFromIDATs <- function(...) {

  targets <- getSamps() 
  rgSet <- read.metharray.exp(".", targets=targets)
  colnames(rgSet) <- rgSet$subject
  grSet <- sesamize(rgSet, ...)
  metadata(grSet)$control_metrics <- t(control_metrics(rgSet)) 
  metadata(grSet)$control_flagged <- t(flag_control_failures(rgSet)) 
  return(grSet)

}
