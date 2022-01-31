#' another oldie but goodie
#' 
#' plot SNPs living in metadata(x)$SNPs, provided they match sampleNames(x)
#' 
#' @param x       a grSet with SNPs in its metadata(), or a SNP matrix 
#' @param rotate  rotate the subjects onto the side? (FALSE)
#' @param ...     other arguments passed on to Heatmap
#' 
#' @import ComplexHeatmap
#' 
#' @export
plotSNPs <- function(x, rotate=FALSE, ...) { 
  if (is(x, "GenomicRatioSet") | is(x, "RGChannelSet")) { 
    SNPs <- as.matrix(metadata(x)$SNPs) # handle DelayedArray SNPs
    if (is.null(SNPs) | !all(colnames(x) %in% colnames(SNPs))) {
      stop("Your SNPs don't match your samples. Aborting.")
    }
    SNPs <- SNPs[, colnames(x)]
  } else { 
    SNPs <- x
  }
  if (rotate) SNPs <- t(SNPs)
  Heatmap(SNPs, name="BAF", ...)
}
