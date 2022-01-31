#' compute the fraction of sample-specific (not masked) NA probes
#'
#' @param   grSet   a GenomicRatioSet 
#' 
#' @return          the fraction NA (of probes that aren't 100% NA) per column
#' 
#' @export
NAfrac <- function(grSet) {

  na_mat <- is.na(getBeta(grSet))
  rowfrac <- rowSums(na_mat) / ncol(na_mat)
  masked_loci <- sum(rowfrac == 1)
  message(masked_loci, " are masked in all samples; omitting from NA tallies.")
  na_frac <- (colSums(na_mat) - masked_loci) / (nrow(na_mat) - masked_loci)
  worstsample <- colnames(na_mat)[which.max(na_frac)]
  message("Failed probes: median of ", round(median(na_frac * 100), 1), "%, ", 
          "maximum of ", round(max(na_frac * 100), 1), "% ",
          "(in ", worstsample, ")")
  return(na_frac) 

}
