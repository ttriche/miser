# Exported classes -------------------------------------------------------------
setClass("MethArraySet", 
         representation(platform = "character"),
         contains = "RangedSummarizedExperiment")

# Exported functions -----------------------------------------------------------
MethArraySet <- function(Beta = NULL, gr = NULL, platform = NULL, ...) {

  if (is.null(Beta) && is.null(gr)) {
    stop("Need either 'gr' or 'Beta' or both to construct a MethArraySet")
  }
  
  if (is.null(platform)) {
    if (!is.null(Beta)) platform <- inferPlatformFromProbeIDs(rownames(Beta))
    if (!is.null(gr)) platform <- inferPlatformFromProbeIDs(names(gr))
  }
  
  if (is.null(gr)) {
    
    df <- sesameDataGet(paste0(platform, ".address"))
    g <- setdiff(names(df), c("ordering", "controls"))
    gr <- df[[g]]
    
    mappable <- names(gr)
    mappable <- intersect(names(gr), rownames(Beta))
    Beta <- Beta[mappable, ]
    gr <- gr[mappable]
    genome(gr) <- g

  } else {

    stopifnot(all(rownames(Beta) %in% names(gr)))
    gr <- gr[rownames(Beta)]

  }
  
  if (is.null(Beta)) Beta <- matrix(NA_real_, nrow=length(gr))
  assays <- SimpleList(Beta = Beta)
  
  new("MethArraySet", 
      platform = platform, 
      SummarizedExperiment(assays = assays, rowRanges=gr, ...))

}


# Exported methods -------------------------------------------------------------
setMethod("show", signature(object = "MethArraySet"), function(object) {
    callNextMethod()
    .show.platform(object)
    .show.genome(object)
})


setMethod("getBeta", signature(object = "MethArraySet"), function(object) {
    return(assay(object, "Beta"))
})


setMethod("getM", signature(object = "MethArraySet"), function(object) {
    return(logit(assay(object, "Beta")))
})


setMethod("annotation",
    signature(object = "MethArraySet"),
    function(object) object@platform
)


.show.platform <- function(object, indent = "  ") {
  cat(sprintf("%splatform: %s\n", indent, object@platform))
}


.show.genome <- function(object, indent = "  ") {
  g <- unique(genome(object))
  if (length(g) == 1) {
    if (is.na(g)) g <- "unspecified"
    cat(sprintf("%sgenome: %s\n", indent, g))
  } else {
    sapply(g, function(gi) cat(sprintf("%s%s: %s\n", indent, gi)))
  }
}
