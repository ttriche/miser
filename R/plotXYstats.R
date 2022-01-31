#' simple function for plotting XY stats (split by inferred sex if found)
#' 
#' @param grSet a GenomicRatioSet
#' @param ...   additional parameters to pass to ComplexHeatmap::Heatmap
#' 
#' @return a Heatmap
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotXYstats <- function(grSet, ...) {

  if ("XYstats" %in% names(metadata(grSet))) {
    XY <- metadata(grSet)$XYstats
  } else {
    XY <- XYstats(grSet)
  }
  if ("inferred_sex" %in% names(colData(grSet))) {
    Heatmap(XY, name="frac", col=jet, column_split=grSet$inferred_sex, ...)
  } else {
    Heatmap(XY, name="frac", col=jet, ...)
  } 

}
