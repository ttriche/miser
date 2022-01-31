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

  XYpct <- XYstats(grSet) * 100
  if ("inferred_sex" %in% names(colData(grSet))) {
    Heatmap(XYpct, name="%", col=jet, column_split=grSet$inferred_sex, ...)
  } else {
    Heatmap(XYpct, name="%", col=jet, ...)
  } 

}
