#' plot heatmap, remind user to check SNPs 
#' 
#' @param   x   a grSet or matrix of beta values
#' @param   k   how many features to keep (max)
#' 
#' @import      ComplexHeatmap 
#' 
#' @export
firstPass <- function(x, k=500, ...) {
  
  if (is(x, "GenomicRatioSet")) {
    k <- intersect(seqlevels(x), paste0("chr", 1:22))
    x <- getBeta(keepSeqlevels(x, k, pruning.mode="coarse"))
  } 
  toPlot <- byExtremality(x, k=k)
  message("Don't forget to plotSNPs() as well!") 
  Heatmap(as(toPlot, "matrix"), col=jet, name="Methylation", 
          clustering_distance_columns="manhattan",
          clustering_method_columns="ward.D2",
          clustering_distance_rows="manhattan",
          clustering_method_rows="ward.D2",
          ...)

}
