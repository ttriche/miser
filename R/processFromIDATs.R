#' read in a directory's worth of IDATs, flag controls, and sesamize them all
#' 
#' Note: this function simply reads all IDATs in the current directory, 
#'       runs QC, stuffs that into metadata, and sesamizes the result. 
#'
#' @param ...     options to pass to sesame::sesamize
#' @param addgeo  optional: try to annotate from GEO? (FALSE) 
#'
#' @return a GenomicRatioSet (or an rgSet in case of failure)
#' 
#' @export 
processFromIDATs <- function(..., addgeo=FALSE) {

  targets <- getSamps() 
  rgSet <- read.metharray.exp(".", targets=targets)

  # QC (document how the SNPs are done...) 
  colnames(rgSet) <- rgSet$subject
  metadata(rgSet)$SNPs <- getSnpBeta(rgSet)
  metadata(rgSet)$control_metrics <- t(control_metrics(rgSet)) 
  metadata(rgSet)$control_flagged <- t(flag_control_failures(rgSet)) 

  # annotate?
  if (addgeo) { 
    message("Attempting to pull metadata...")
    res <- try(addCharacteristics(rgSet))
    if (inherits(res, "try-error")) {
      message("Failed, returning unmodified results.")
    } else {
      message("Success!")
      rgSet <- res
    }
  }

  # create a masked GenomicRatioSet
  grSet <- try(sesamize(rgSet, ...))
  if (inherits(grSet, "try-error")) {
    message("Failed to create GenomicRatioSet! Returning raw RGChannelSet.")
    return(rgSet)
  } else { 
    metadata(grSet) <- metadata(rgSet)
    colData(grSet)$inferred_sex <- getSex(grSet)
    # sesame's version of this is better but for now just use minfi's 
    return(grSet)
  }
  # eventually, it would be better to create a grSet-like object that wraps 
  # a FileSet but links to appropriate metadata, e.g. a GenomicFileSet 

}
