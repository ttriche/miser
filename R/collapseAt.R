#' Collapse a SummarizedExperiment-like object over a GRanges
#'
#' @param x         anything descended from a SummarizedExperiment
#' @param y         regions to summarize over, as a GRanges
#' @param how       medians or means? (medians)
#' @param minprobes minimum probes to summarize over (1)
#' @param imp       impute any remaining NAs? (FALSE) 
#' @param obj       return an object resembling `x`? (FALSE) 
#'
#' @return          an object with same colData but with new rowRanges & assays
#' 
#' @details
#' No imputation is performed and only unique ranges are summarized.
#' If you want to disjoin ranges, do it before running this function. 
#' The function just takes whatever is in assays(x)[[1]] and summarizes that.
#' 
#' @import          minfi  
#' @import          impute 
#' @import          GenomicRanges 
#'
#' @export
collapseAt <- function(x, y, how=c("medians", "means"), minprobes=2, imp=TRUE, obj=FALSE) {

  how <- match.arg(how) 
  fn <- switch(how,
               medians=colMedians,
               means=colMeans)
  if (class(y) == "dmrcate.output") y <- extractRanges(y)
  stopifnot(is(y, "GenomicRanges")) ## must be a GRanges
  stopifnot(is(x, "SummarizedExperiment"))

  y <- unique(y)
  names(y) <- as.character(y)

  if (is(x, "GenomicRatioSet")) x <- x[ grep("^cg", rownames(x)), ] 

  ## find and index the runs
  yy <- subsetByOverlaps(y, x)
  xx <- subsetByOverlaps(x, yy)
  xx <- xx[rowSums(is.na(assay(xx, 1))) < (ncol(x)/2), ]
  yy <- subsetByOverlaps(yy, xx)
  ol <- findOverlaps(xx, yy)

  byDMR <- as(lapply(split(as(ol, "data.frame"), subjectHits(ol)),
                     function(z) z[, "queryHits"]), "List")
  names(byDMR) <- names(yy)
  if (minprobes > 1) byDMR <- byDMR[sapply(byDMR, length) >= minprobes]
    
  summarize <- 
    switch(how,
           medians=function(w, b) colMedians(b[w,,drop=F], na.rm=TRUE),
           means=function(w, b) colMeans(b[w,,drop=F], na.rm=TRUE))
  res <- do.call(rbind, lapply(byDMR, summarize, b=as.matrix(getBeta(xx))))
  NAs <- length(which(is.na(res)))
  total <- prod(dim(res))
  pctNA <- (NAs/total) * 100

  message(NAs, "/", total, " (", round(pctNA, 3), "%) NAs after collapsing.")
  if (imp) {
    message("Imputing... ")
    res <- impute.knn(res)$data
  }

  colnames(res) <- colnames(xx)

  if (!obj) {
    return(res)
  } else { 
    xx <- x[seq_len(nrow(res)),] 
    rownames(xx) <- rownames(res)
    rr <- as(rownames(res), "GRanges")
    rr$score <- countOverlaps(rr, x)
    rowRanges(xx) <- rr
    assays(xx) <- assays(xx)[1]
    assays(xx, withDimnames=FALSE)[[1]] <- res
    rownames(xx) <- rownames(res)
    return(xx)
  } 
}
