#' compute assorted sex chromosome related statistics (XNA, YNA, XBeta)
#'
#' @param   grSet   a GenomicRatioSet 
#' 
#' @return          a matrix of XNA, YNA, XBeta for each column
#' 
#' @export
XYstats <- function(grSet) {

  XNA <- is.na(getBeta(subset(grSet, seqnames %in% c("chrX", "X"))))
  rowfrac <- rowSums(XNA) / ncol(XNA)
  masked_X <- sum(rowfrac == 1)
  XNAfrac <- (colSums(XNA) - masked_X) / (nrow(XNA) - masked_X)

  YNA <- is.na(getBeta(subset(grSet, seqnames %in% c("chrY", "Y"))))
  yrowfrac <- rowSums(YNA) / ncol(YNA)
  masked_Y <- sum(yrowfrac == 1)
  YNAfrac <- rep(1, ncol(YNA))
  if (!all(yrowfrac == 1)) {
    YNAfrac <- (colSums(YNA)-masked_Y)/(nrow(YNA)-masked_Y)
  }

  if (!"IslandStatus" %in% names(rowData(grSet))) {
    rowData(grSet)$IslandStatus <- minfi::getIslandStatus(grSet)
  }
  whichX <- rownames(subset(grSet, seqnames %in% c("chrX", "X") & 
                                   rowData(grSet)$IslandStatus != "OpenSea"))
  XBeta <- colMedians(getBeta(grSet[whichX, ]), na.rm=TRUE)
  
  return(rbind(XNAfrac=XNAfrac, YNAfrac=YNAfrac, XBeta=XBeta))

}
