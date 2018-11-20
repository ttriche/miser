#' punishingly simple VAF clustering for clonality estimation
#' 
#' Suppose we have read coverage over a base with a detected variant.
#' We can model the reads mapping to variant and reference bases as 
#' 
#'   p ~ Beta(variant_reads, reference_reads)
#' 
#' Now, if we take logits, with nu = logit(p), then
#' 
#'   nu ~ Normal( log(p/(1-p)), log10(total reads) )
#' 
#' As a bonus, we can model `nu` as roughly multivariate normal.
#' Furthermore, we can use the latter term as observation weights.
#' Given a vector or matrix of VAFs, this function clusters them
#' using dbscan, with log10(read_depth) as the weight for each VAF.
#' 
#' @param   alt   read depth for the variant allele 
#' @param   ref   read depth for the reference allele 
#' @param   eps   epsilon to use for DBSCAN clustering 
#' @param   ...   additional arguments to pass on to dbscan 
#' 
#' @return        an object of class `dbscan_fast`
#' 
#' @import dbscan
#'
#' @export
vafClust <- function(alt, ref, depth=NULL, ...) { 
 
  stopifnot(identical(dim(alt), dim(ref)))
  total <- alt + ref
  VAF <- alt / total 
  dbscan::dbscan(x=flogit(VAF), eps=eps, weights=log10(total), ...)

}
