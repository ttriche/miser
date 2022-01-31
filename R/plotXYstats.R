#' simple function for plotting XY stats (split by inferred sex if found)
#' 
#' @param grSet   a GenomicRatioSet
#' @param column  optional column name with sex information
#' @param ...     additional parameters to pass to ComplexHeatmap::Heatmap
#' 
#' @return a Heatmap
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotXYstats <- function(grSet, column=NULL, ...) {

  if ("XYstats" %in% names(metadata(grSet))) {
    XY <- metadata(grSet)$XYstats
  } else {
    XY <- XYstats(grSet)
  }
  if (!is.null(column) & column %in% names(colData(grSet))) {
    message("Plotting based on `", column, "`...")
    Heatmap(XY, name="frac", col=jet, column_split=colData(grSet)[,column], ...)
  } else if ("sex" %in% names(colData(grSet))) {
    message("Plotting based on `sex`...")
    Heatmap(XY, name="frac", col=jet, column_split=grSet$sex, ...)
  } else if ("inferred_sex" %in% names(colData(grSet))) {
    message("Plotting based on `inferred_sex`...")
    Heatmap(XY, name="frac", col=jet, column_split=grSet$inferred_sex, ...)
  } else {
    Heatmap(XY, name="frac", col=jet, ...)
  } 

}
