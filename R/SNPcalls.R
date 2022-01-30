#' fit a 1-to-3-component mixture model to flogit(SNPbetas) for each sample
#' 
#' @param SNPs      a GenomicRatioSet or a matrix of SNP beta values 
#'
#' @return          0/1/2 calls for each SNP in each sample
#' 
#' @details         Mclust(miser::flogit(x), G=1:3) is run for each sample. 
#' 
#' @import mclust
#' @import minfi
#' 
#' @export
SNPcalls <- function(SNPs) {

  if (is(SNPs, "GenomicRatioSet") | is(SNPs, "RGChannelSet")) {
    if (!"SNPs" %in% names(metadata(SNPs))) stop("No `SNPs` in metadata.")
    SNPs <- metadata(SNPs)$SNPs
  }

  # Could be parallelized if truly absurd sample sizes are being fit. 
  apply(SNPs, 2, .SNPcall)

}


# helper for heavy lifting
.SNPcall <- function(SNPs) {

  SNPfit <- Mclust(as.numeric(flogit(SNPs)), G=1:3)
  calls <- predict(SNPfit, as.numeric(flogit(SNPs)))$classification - 1
  names(calls) <- names(SNPs)
  return(calls) 

}
