#' read in a directory's worth of IDATs, flag controls, sesamize, then dump
#' 
#' for ease of factorization, inferred sex and flagged controls are written out
#' as a control/QC metrics matrix, along with the fraction of NAs among probes
#' that would not typically be masked by sesame, and the median total intensity
#' 
#' @param stub    mandatory prefix for files (stub_betas.csv, stub_qc.csv)
#' @param path    where to save the dumped CSV files (".")
#' @param ...     options to pass to sesamize
#'
#' @return a GenomicRatioSet (or an rgSet in case of failure)
#' 
#' @export 
dumpFromIDATs <- function(stub, path=".", ...) {

  grSet <- try(processFromIDATs(..., addgeo=addgeo))
  if (!is(grSet, "GenomicRatioSet")) stop("Failed to sesamize to a grSet.")
  dumpQCfiles(grSet, stub=stub, path=path, betas=TRUE)
  return(grSet)

}
