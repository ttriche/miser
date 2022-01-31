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

  message("Attempting to sesamize ", ncol(rgSet), " samples...") 
  grSet <- try(sesame::sesamize(rgSet, ...))
  if (inherits(grSet, "try-error")) {
    message("Failed to create GenomicRatioSet! Returning raw RGChannelSet.")
    return(rgSet)
  } else { 
    message("Done.")
  }

  metadata(grSet) <- metadata(rgSet)
  if (!"SNPs" %in% names(metadata(grSet))) {
    message("Adding SNP intensities to metadata...")
    metadata(grSet)$SNPs <- getSnpBeta(rgSet) 
  } 
  if (!"control_flagged" %in% names(metadata(grSet))) {
    message("Flagging control probe failures in metadata...")
    metadata(grSet)$control_flagged <- t(flag_control_failures(rgSet))
  }
  # sesame's version is better,
  # but minfi's is more convenient
  infsex <- try(minfi::getSex(grSet))
  if (!inherits(infsex, "try-error")) {
    colData(grSet)$inferred_sex <- infsex$predictedSex
  } else { 
    warning("Could not automatically infer sex from CN.")
  }
  rowData(grSet)$IslandStatus <- minfi::getIslandStatus(grSet)
  metadata(grSet)$XYstats <- XYstats(grSet) # using the above bits...
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
