#' call and plot SNPs from a matrix of SNP betas 
#' 
#' plot SNPs living in metadata(x)$SNPs, provided they match sampleNames(x)
#' 
#' @param x       a grSet with SNPs in its metadata()
#' @param rotate  rotate the subjects onto the side? (FALSE)
#' @param ...     other arguments passed on to Heatmap
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotSNPcalls <- function(x, rotate=FALSE, ...) { 

  SNP <- colorRamp2(seq(0, 2), c("#00007F", "yellow", "#7F0000"))
  SNPs <- as.matrix(metadata(x)$SNPs) # handle DelayedArray SNPs
  if (is.null(SNPs) | !all(colnames(x) %in% colnames(SNPs))) {
    stop("Your SNPs don't match your samples. Aborting.")
  }
  SNPcalls <- SNPcalls(SNPs[, colnames(x)])
  if (rotate) SNPcalls <- t(SNPcalls)
  Heatmap(SNPcalls, col=SNP, name="Alleles", ...)

}
