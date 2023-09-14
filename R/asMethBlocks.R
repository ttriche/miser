#' Collapse a GenomicRatioSet or similar over methylation blocks (per Kaplan)
#'
#' @param x         a GenomicRatioSet or something like it
#' @param genome    genome to use for methylation block coordinates (hg19)
#' @param proper    restrict to methylation blocks defined in Kaplan? (FALSE)
#' @param stable    restrict to methylation blocks shared by hg19&hg38? (FALSE)
#'
#' @return          an object with same colData but with new rowRanges & assays
#' 
#' @details
#' By default, data(methBlocks) is used for summarization, singleton probes are
#' included, methylation block rownames are returned against hg19, and neither 
#' singletons (!proper) nor hg19-only blocks (!stable) are filtered out.
#' 
#' @import          minfi  
#' @import          GenomicRanges 
#'
#' @export
asMethBlocks <- function(x, genome="hg19", proper=FALSE, stable=FALSE) {

  stopifnot(is(x, "SummarizedExperiment"))
  
  data("methBlocks", package="miser") 
  keepCpGs <- rownames(x)[rowSums(is.na(getBeta(x))) < ncol(x)]
  methBlocks <- subset(methBlocks, rownames(methBlocks) %in% keepCpGs)
  prop <- paste0("in", toupper(substr(genome, 1, 1)), substr(genome, 2, 4))
  if (proper) methBlocks <- subset(methBlocks, methBlocks[[prop]])
  if (stable) methBlocks <- subset(methBlocks, inHg19 & inHg38)
  methBlocks <- subset(methBlocks, !is.na(methBlocks[[genome]]))
  message("Computing per-methylation-block averages...")
  blockBetas <- do.call(rbind, 
                        by(getBeta(x[rownames(methBlocks),]), 
                           INDICES=methBlocks[[genome]],
                           FUN=colMeans, na.rm=TRUE))
  message("Done. Reconstructing ", class(x), "...")
  amb <- x[seq_len(nrow(blockBetas)), ]
  rownames(amb) <- rownames(blockBetas)
  assay(amb, "Beta") <- blockBetas
  rowRanges(amb) <- as(rownames(blockBetas), "GRanges")
  genome(amb) <- unique(genome(amb))
  return(amb)

}
