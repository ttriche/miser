#' idiotically simple masking and imputation
#' 
#' @param x   a GenomicRatioSet
#' @param na  probes with this proportion NA (or higher) are dropped (0.5)
#' @param sqz a squeeze factor to regularize imputation (1 - 1e-6 by default) 
#' 
#' @return    a GenomicRatioSet with masked probes dropped and NAs imputed 
#' 
#' @import    impute
#' @import    minfi 
#' 
#' @export
fixNAs <- function(x, na=0.5, sqz=(1 - 1e-6)) { 
  stopifnot(is(x, "GenomicRatioSet")) 
  stopifnot("mask" %in% names(rowData(x)))
  naFrac <- rowSums(is.na(getBeta(x)))
  rowData(x)$mask <- rowData(x)$mask | naFrac >= na
  x <- subset(x, !rowData(x)$mask)
  if (any(is.na(getBeta(x)))) {
    assays(x)$Beta <- inv.logit(impute.knn(logit(getBeta(x) * sqz))$data)/sqz
  }
  return(x) 
}
