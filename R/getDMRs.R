#' use DMRcate to get DMRs (with less hassle) 
#' 
#' This is a convenience wrapper that leaves out much of DMRcate's flexibility.
#' If you need that flexibility, you are well advised to use DMRcate directly. 
#' 
#' @param x               a GenomicRatioSet 
#' @param design          a design matrix
#' @param dropXY          drop sex chromosomes? (TRUE, usually a good idea)
#' @param impute          impute NAs (FALSE; this will hang on HDF5 if TRUE)
#' @param coef            which column to fit (default is 2)
#' @param betacutoff      DMRs must have at least this maxbetaFC (0.1)
#' @param ...             other arguments to pass to DMRcate::cpg.annotate
#' 
#' @return                tidied output from DMRcate::dmrcate
#' 
#' @import DMRcate
#'
#' @export
getDMRs <- function(x, design, dropXY=TRUE, impute=FALSE, coef=2, fdr=.05, 
                    betacutoff=.1, DMLs=FALSE, parallel=FALSE, ...) {

  message("Note: this is a simplified, differential-only version of DMRcate.")
  message("DMRcate is capable of much more involved analyses. Read its manual.")
  stopifnot(!is.null(design))

  message("Dropping non-CpG probes...")
  x <- subset(x, substr(rownames(x), 1, 2) == "cg")
  if (dropXY) {
    message("Dropping sex chromosomes if present (dropXY was set to TRUE)...")
    x <- keepSeqlevels(x, paste0("chr", 1:22), pruning.mode="coarse")
  }
  
  # avoid hanging with HDF5 SEs?
  if (impute) {
    message("Imputing (this may take a VERY LONG TIME with HDF5 backing)...")
    x <- fixNAs(x) # the alternative is to fail...
    message("Done imputing.")
  }
  if ("mask" %in% names(rowData(x))) x <- subset(x, !rowData(x)$mask)
  
  message("Annotating individual CpGs...")
  annot <- cpg.annotate("array", x, design=design, coef=2, fdr=fdr, ...)
 
  message("Demarcating significant regions...")
  res <- extractRanges(dmrcate(annot, C=2, betacutoff=betacutoff), 
                       genome=unique(genome(x)))
  message("Found ", length(res), " DMRs.")
  return(res)

}
