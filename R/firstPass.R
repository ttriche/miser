#' plot heatmap, remind user to check SNPs 
#' 
#' @param   x   a grSet
#' @param   k   how many features to keep (max)
#' 
#' @return      a Heatmap object from ComplexHeatmap 
#' 
#' @import      ComplexHeatmap 
#' 
#' @export
firstPass <- function(x, k=500, ...) {

  toPlot <- byExtremality(getBeta(keepSeqlevels(x, paste0("chr", 1:22), 
                                                pruning.mode="coarse")), k=k)
  message("Don't forget to plotSNPs() as well!") 
  Heatmap(toPlot, col=jet, name="Methylation", 
          clustering_distance_columns="manhattan",
          clustering_method_columns="ward.D2",
          clustering_distance_rows="manhattan",
          clustering_method_rows="ward.D2",
          ...)

}



