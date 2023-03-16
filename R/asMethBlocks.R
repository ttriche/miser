#' Collapse a GenomicRatioSet or similar over methylation blocks (per Kaplan)
#'
#' @param x         usually, a GenomicRatioSet or something like it
#' @param y         regions to summarize over (methBlocks.EPIC.hg19)
#' @param useCpH    include CpH probes? (FALSE, because it's a bad idea)
#' @param minprobes minimum probes to summarize over (1)
#' @param obj       keep the object the same (e.g. summarizedExperiment)? (TRUE)
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
asMethBlocks <- function(x, y=NULL, useCpH=FALSE, minprobes=1, obj=TRUE) {

  require(matrixStats)
  if (is.null(y)) {
    mb <- paste0("methBlocks.EPIC.", unique(genome(x)))
    data(mb, package="miser")
    y <- get(mb)
  }
  stopifnot(is(y, "GenomicRanges")) ## must be a GRanges
  stopifnot(is(x, "SummarizedExperiment"))
  if (is(x, "GenomicRatioSet") & !useCpH) x <- x[ grep("^cg", rownames(x)), ] 
  collapseAt(x=x, y=y, minprobes=minprobes, imp=FALSE, obj=obj)

}
