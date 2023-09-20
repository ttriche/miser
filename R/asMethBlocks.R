#' Collapse a GenomicRatioSet or similar over methylation blocks (per Kaplan)
#' 
#' By default, data(methBlocks) is used for summarization, singleton probes are
#' included, methylation block rownames are returned against hg19, and neither 
#' singletons (!proper) nor unstable mappings (!stable) are filtered out.
#'
#' @param x         a GenomicRatioSet or something like it
#' @param g         genome to use for methylation block coordinates (hg19)
#' @param proper    restrict to methylation blocks defined in Kaplan? (FALSE)
#' @param stable    restrict to sites/blocks shared across hg19 & hg38? (FALSE)
#'
#' @return          an object with same colData but with new rowRanges & assays
#' 
#' @details
#' `stable` is a weird flag and you probably should not use it. If and only if 
#' a probe maps to the same region in hg19 and hg38, it is deemed stable. That
#' region can be one base pair wide, or it can be a thousand base pairs wide. 
#' The point is that this set of mappings is very small (about 4000 probes). 
#' 
#' @seealso         switchMethBlocksGenome
#' 
#' @import          minfi  
#' @import          BiocParallel 
#' @import          GenomicRanges 
#'
#' @export
asMethBlocks <- function(x, g=c("hg19","hg38"), proper=FALSE, stable=FALSE) {

  g <- match.arg(g)
  g0 <- setdiff(c("hg19","hg38"), g)
  stopifnot(is(x, "SummarizedExperiment"))
  if (exists("methBlocks")) {
    message("You seem to have defined `methBlocks` already. Using that.")
  } else { 
    data("methBlocks", package="miser") 
  }

  message("Checking for masked (NA) probes...")
  N <- ncol(x)
  keepCpGs <- rownames(x)[rowSums(is.na(getBeta(x))) < N]
  M <- length(keepCpGs)
  dropCpGs <- setdiff(rownames(x), keepCpGs)
  M_NA <- length(dropCpGs)
  message("Retained ", M, " CpGs, dropped ", M_NA, " with ", N, "/", N, " NAs.")

  message("Mapping probes to methylation blocks in ", g, "...")
  methBlocks <- methBlocks[keepCpGs, ]
  methBlocks <- subset(methBlocks, !is.na(methBlocks[[g]]))
  prop <- paste0("in", toupper(substr(g, 1, 1)), substr(g, 2, 4), "mb")
  if (proper) {
    stopifnot(prop %in% names(methBlocks))
    methBlocks <- subset(methBlocks, methBlocks[[prop]])
    keepCpGs <- intersect(keepCpGs, rownames(methBlocks))
    M <- length(keepCpGs)
    message("Retaining (only) ", M, " probes in defined blocks.")
    message("If this is not what you expected, set `proper` to FALSE.")
  } else { 
    message("Retaining singleton probes located outside of defined blocks.")
  }
  g0 <- setdiff(c("hg19","hg38"), g)
  if (stable) {
    methBlocks <- subset(methBlocks, stable) # see details
    keepCpGs <- intersect(keepCpGs, rownames(methBlocks))
    M <- length(keepCpGs)
    message("Retaining (only) ", M, " stably mapped probes.")
    message("If this is not what you expected, set `stable` to FALSE.")
  } else {
    keepCpGs <- intersect(keepCpGs, rownames(methBlocks))
    M <- length(keepCpGs)
    message("Retained ", M, " probes mapped to ", g, " (and most to ", g0, ").")
  }

  # preparation for aggregation
  tbl <- table(methBlocks[[g]])
  blocks <- names(tbl)[tbl > 1] # blocks with more than one probe
  byBlock <- subset(methBlocks, methBlocks[[g]] %in% blocks)[, g, drop=FALSE]

  # this will become rowRanges(amb)
  mbgr <- as(methBlocks[[g]], "GRanges")
  genome(mbgr) <- g 
  mbgr$probes <- countOverlaps(mbgr, granges(x)[keepCpGs])
  mbgr$singleton <- width(mbgr) == 1
  mbgr <- unique(mbgr)
  names(mbgr) <- as.character(mbgr)

  message("Computing per-block average methylation across ", M, " probes...")
  placeholder <- rownames(subset(methBlocks, !duplicated(methBlocks[[g]])))
  blockBetas <- getBeta(x[placeholder, ])
  rownames(blockBetas) <- names(mbgr) 
  # only recompute regions with more than one probe
  byBlock[[g]] <- factor(byBlock[[g]])
  toSquash <- split.data.frame(getBeta(x[rownames(byBlock), ]), byBlock[[g]])
  res <- do.call(rbind, lapply(toSquash, colMeans, na.rm=TRUE)) # bplapply fails
  blockBetas[rownames(res), ] <- res
 
  message("Reconstructing ", class(x), "...")
  amb <- x[seq_len(nrow(blockBetas)), ]
  rownames(amb) <- rownames(blockBetas)
  assay(amb, "Beta") <- blockBetas
  rowRanges(amb) <- mbgr
  
  message("Done.")
  return(amb)

}


# helper function
sharedMethBlockPairs <- function(mb, g=c("hg19","hg38")) {

  g <- match.arg(g)
  g0 <- "hg38"
  mb <- subset(mb, !is.na(mb[[g]]) & !is.na(mb[[g0]]))
  mb$inBlocks <- mb[["inHg19mb"]] & mb[["inHg38mb"]]
  mb$pairing <- paste0(as.character(mb[[g]]), "|", as.character(mb[[g0]]))
  dupes <- mb$pairing[which(duplicated(mb$pairing))]
  mb$dupe <- mb$pairing %in% dupes
  mb$hits <- 1
  mb$hits <- rowsum(mb$hits, mb$pairing)[mb$pairing, 1]
  mb <- subset(mb, !duplicated(mb$pairing))

  message("Pairing ", g, " and ", g0, " methylation block ranges...")
  paired <- Pairs(as(mb[[g]], "GRanges"), as(mb[[g0]], "GRanges"))
  names(paired) <- as.character(mb[[g]])
  mcols(paired)$inBlocks <- mb$inBlocks
  dupes <- names(paired)[duplicated(names(paired))]
  mcols(paired)$splitMap <- names(paired) %in% dupes
  ndups <- sum(names(paired) %in% dupes)
  nuniq <- length(paired) - ndups
  message("Found ", ndups, " split mappings and ", 
          nuniq, " unique mappings from ", g, " to ", g0, ".")
  return(paired)

}

# helper function
splitSharedMethBlockPairs <- function(mb, g=c("hg19","hg38")) {

  paired <- sharedMethBlockPairs(mb, g)
  splt <- subset(paired, mcols(paired)$splitMap)
  ssplt <- split(as.character(second(splt)), as.character(first(splt)))

}
