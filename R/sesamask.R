#' add a mask to rowData(grSet) in a platform-appropriate way
#' 
#' since GenomicRatioSets have their platform annotated, we just use that 
#' 
#' @param x   a suitable object (SummarizedExperiment-derived)
#' 
#' @return    the same object but with $mask in rowData(x) filled out
#' 
#' @import sesame
#' @import sesameData
#' 
#' @export
sesamask <- function(x) {
  if (is.null(annotation(x)["array"])) stop("Cannot determine mask to use.") 
  platform <- switch(annotation(x)['array'],
                     IlluminaHumanMethylationEPIC="EPIC.probeInfo",
                     IlluminaHumanMethylation450k="HM450.probeInfo",
                     IlluminaHumanMethylation27k="HM27.probeInfo")
  rowData(x)$mask <- rownames(x) %in% sesameDataGet(platform)$mask
  return(x)
}
