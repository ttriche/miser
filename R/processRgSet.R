#' sesamize a result with sensible QC and passing BPPARAM, etc.
#' 
#' @param x       a SummarizedExperiment-derived object with SNPs and controls
#' @param justDump  just dump the sesamized result? (TRUE) 
#' @param addgeo    optional: try to annotate from GEO, if not already done? (F)
#' @param ...       options to pass to sesamize
#'
#' @return          a GenomicRatioSet (or an res if failure)
#'
#' @import sesame 
#' @import minfi
#'
#' @export 
processResult <- function(x, addgeo=FALSE, justResult=TRUE, ...) {

  message("Attempting to sesamize ", ncol(x), " samples...") 
  res <- try(sesamize(x, ...)) # now part of miser
  if (inherits(res, "try-error")) {
    message("Failed to create GenomicRatioSet! Returning raw RGChannelSet.")
    return(x)
  } else { 
    message("Done.")
  }
  metadata(res) <- metadata(x)

  if (!justDump) {
    # {{{ QC and cleanup as needed
    if (!"SNPs" %in% names(metadata(res))) {
      message("Adding SNP intensities to metadata...")
      metadata(res)$SNPs <- getSnpBeta(x) 
    } 
    if (!"control_flagged" %in% names(metadata(res))) {
      message("Flagging control probe failux in metadata...")
      metadata(res)$control_flagged <- t(flag_control_failures(x))
    }
    # sesame's version is better,
    # but minfi's is more convenient
    infsex <- try(minfi::getSex(res))
    if (!inherits(infsex, "try-error")) {
      colData(res)$inferred_sex <- infsex$predictedSex
    } else { 
      warning("Could not automatically infer sex from CN.")
    }
    rowData(res)$IslandStatus <- minfi::getIslandStatus(res)
    metadata(res)$XYstats <- XYstats(res) # using the above bits...
    rowData(res)$NAfrac <- rowSums(is.na(getBeta(res)))/ncol(res)
    colData(res)$NAfrac <- NAfrac(res)
    # }}}
    if (addgeo) { 
      message("Attempting to pull metadata...")
      # {{{
      annotated <- try(addCharacteristics(res))
      if (inherits(annotated, "try-error")) {
        message("GEO annotation failed, returning unmodified results.")
      } else {
        message("Success! ", ncol(annotated), " samples annotated from GEO.")
        res <- annotated
      }
      # }}}
    }
  }

  return(res)

}


# backwards compatibility
processRgSet <- processResult
