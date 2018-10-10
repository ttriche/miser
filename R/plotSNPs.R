#' another oldie but goodie
#' 
#' plot SNPs living in metadata(x)$SNPs, provided they match sampleNames(x)
#' 
#' @param x   a grSet with SNPs in its metadata()
#' 
#' @import ComplexHeatmap
#' 
#' @export
plotSNPs <- function(x) { 
  SNPs <- metadata(x)$SNPs
  if (is.null(SNPs) | !all(colnames(x) %in% colnames(SNPs))) {
    stop("Your SNPs don't match your samples. Aborting.")
  }
  Heatmap(SNPs[, colnames(x)], name="BAF")
}
