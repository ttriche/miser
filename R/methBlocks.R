#' Collapse a GenomicRatioSet or similar over methylation blocks (per Kaplan)
#'
#' @param x         anything descended from a SummarizedExperiment
#' @param y         regions to summarize over (methBlocks.EPIC.hg19)
#' @param useCpH    include CpH probes? (FALSE, because it's a bad idea)
#' @param minprobes minimum probes to summarize over (1)
#'
#' @return          an object with same colData but with new rowRanges & assays
#' 
#' @details
#' By default, methBlocks.EPIC.hg19 is used for summarization, and CpHs are not.
#' 
#' @import          minfi  
#' @import          GenomicRanges 
#'
#' @export
methBlocks <- function(x, y=NULL, useCpH=FALSE, minprobes=1) {

  require(matrixStats)
  if (is.null(y)) {
    data("methBlocks.EPIC.hg19")
    y <- methBlocks.EPIC.hg19
  }
  stopifnot(is(y, "GenomicRanges")) ## must be a GRanges
  stopifnot(is(x, "SummarizedExperiment"))
  if (is(x, "GenomicRatioSet") & !useCpH) x <- x[ grep("^cg", rownames(x)), ] 
  collapseAt(x=x, y=y, minprobes=minprobes, imp=FALSE)

}
