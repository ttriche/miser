#' sesamize an RGChannelSet with sensible QC and passing BPPARAM, etc.
#' 
#' @param rgSet     an RGChannelSet with SNPs and control_flagged in metadata()
#' @param addgeo    optional: try to annotate from GEO, if not already done? (F)
#' @param ...       options to pass to sesame::sesamize
#'
#' @return          a GenomicRatioSet (or an rgSet if failure)
#'
#' @import sesame 
#' @import minfi
#'
#' @export 
processRgSet <- function(rgSet, addgeo=FALSE, ...) {
  
  grSet <- try(sesame::sesamize(rgSet, ...))
  if (inherits(grSet, "try-error")) {
    message("Failed to create GenomicRatioSet! Returning raw RGChannelSet.")
    return(rgSet)
  }

  metadata(grSet) <- metadata(rgSet)
  # sesame's version of this is better, but for now just use minfi's 
  colData(grSet)$inferred_sex <- minfi::getSex(grSet)$predictedSex
  rowData(grSet)$IslandStatus <- minfi::getIslandStatus(grSet)
  rowData(grSet)$NAfrac <- rowSums(is.na(getBeta(grSet)))/ncol(grSet)
  colData(grSet)$NAfrac <- NAfrac(grSet)
  
  if (addgeo) { 
    message("Attempting to pull metadata...")
    res <- try(addCharacteristics(grSet))
    if (inherits(res, "try-error")) {
      message("GEO annotation failed, returning unmodified results.")
    } else {
      message("Success! ", ncol(res), " samples annotated with GEO metadata.")
      grSet <- res
    }
  }
  
  return(grSet)

}
