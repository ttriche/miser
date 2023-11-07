#' read in a directory's worth of IDATs, QC them, optionally sesamize them all
#' 
#' Note: this function simply reads all IDATs in the current directory, 
#'       runs QC, stuffs in metadata, and (if !justRgSet) sesamizes it. 
#'
#' @param frags     which elements of the filenames are relevant? (1:3)
#' @param targets   a targets dataframe (will be passed to minfi for reading)
#' @param addgeo    optional: try to annotate from GEO? (FALSE) 
#' @param justRgSet optional: dump the result and don't sesamize? (FALSE)
#' @param forSesame optional: avoid loading minfi data structures? (FALSE) 
#' @param ...       options to pass to sesamize
#'
#' @return a SummarizedExperiment-derived object of some sort (see Details)
#'
#' @details
#' minfi is becoming a pain in the ass to deal with due to its dependencies, so
#' this function can emit various SummarizedExperiment-derived results (whether
#' an RGChannelSet, a GenomicRatioSet, a SummarizedExperiment, or a 
#' SingleCellExperiment, depending upon its parameters). The only real guarantee
#' is that, in the absence of an error, a SummarizedExperiment-derived object 
#' will be returned, and it will accommodate the usual generic getWhatever()s. 
#' 
#' @import R.utils
#' @import sesame 
#' @import minfi
#'
#' @export 
#'
processFromIDATs <- function(frags=1:3, targets=NULL, addgeo=FALSE, justRgSet=FALSE, forSesame=FALSE, verbose=FALSE, ...) {

  if (is.null(targets)) {
    if (verbose) message("Cataloging IDATs...")
    targets <- getSamps(frags=frags) 
    message(nrow(targets), " samples found. Reading signals...")
  } 

  if (!"subject" %in% names(targets)) targets[, "subject"] <- rownames(targets)

  if (forSesame) { 
    stop("minfi-free operation is not yet supported (soon!)...")
  } else { 
    res <- read.metharray.exp(".", targets=targets, verbose=TRUE, force=TRUE)
  }

  # QC (document how the SNPs are done...) 
  colnames(res) <- res$subject
  metadata(res)$SNPs <- getSnpBeta(res)
  metadata(res)$control_metrics <- t(control_metrics(res)) 
  metadata(res)$control_flagged <- t(flag_control_failures(res)) 
  message(ncol(res), " samples QC'ed. Proceeding...")

  # annotate?
  if (addgeo) { 
    message("Attempting to pull metadata...")
    res <- try(addCharacteristics(res))
    if (inherits(res, "try-error")) {
      message("GEO annotation failed, returning unmodified results.")
    } else {
      message("Success!")
      res <- res
    }
  }

  # just return the res?
  if (justRgSet) return(res)

  # create GenomicRatioSet?
  # nb. this should be deprecated
  processResult(res, ...) 

}
