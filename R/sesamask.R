#' add a mask to rowData(grSet) in a platform-appropriate way
#' 
#' since GenomicRatioSets have their platform annotated, we just use that 
#' 
#' @param x   a GenomicRatioSet
#' 
#' @return    the same GenomicRatioSet but with $mask in rowData(x) filled out
#' 
#' @import sesame
#' @import sesameData
#' 
#' @export
sesamask <- function(x) {
  stopifnot(is(x, "GenomicRatioSet"))
  platform <- switch(annotation(x)['array'],
                     IlluminaHumanMethylationEPIC="EPIC.probeInfo",
                     IlluminaHumanMethylation450k="HM450.probeInfo",
                     IlluminaHumanMethylation27k="HM27.probeInfo")
  rowData(x)$mask <- rownames(x) %in% sesameDataGet(platform)$mask
  return(x)
}
