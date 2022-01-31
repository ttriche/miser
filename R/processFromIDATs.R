#' read in a directory's worth of IDATs, QC them, optionally sesamize them all
#' 
#' Note: this function simply reads all IDATs in the current directory, 
#'       runs QC, stuffs in metadata, and (if !justRgSet) sesamizes it. 
#'
#' @param frags     which elements of the filenames are relevant? (1:3)
#' @param addgeo    optional: try to annotate from GEO? (FALSE) 
#' @param justRgSet optional: dump the rgSet and don't sesamize? (FALSE)
#' @param ...       options to pass to sesame::sesamize
#'
#' @return a GenomicRatioSet (or an rgSet in case of failure)
#'
#' @import R.utils
#' @import sesame 
#' @import minfi
#'
#' @export 
processFromIDATs <- function(frags=1:3, addgeo=FALSE, justRgSet=FALSE, ...) {

  message("Cataloging IDATs...")
  targets <- getSamps(frags=frags) 
  message(nrow(targets), " samples found. Reading signals...")
  rgSet <- read.metharray.exp(".", targets=targets, verbose=TRUE)

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
      message("GEO annotation failed, returning unmodified results.")
    } else {
      message("Success!")
      rgSet <- res
    }
  }

  # just return the rgSet?
  if (justRgSet) return(rgSet)

  # create GenomicRatioSet
  processRgSet(rgSet, ...) 

}
