#' simple function for plotting XY stats (split by inferred sex if found)
#' 
#' @param x       a GenomicRatioSet, or an XYstats matrix 
#' @param sex     optional sex labels to use for splitting the plot 
#' @param ...     additional parameters to pass to ComplexHeatmap::Heatmap
#' 
#' @return a Heatmap
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotXYstats <- function(x, sex=NULL, ...) {

  if (!is.null(sex)) message("Splitting plot on provided sex covariate...")
  if (is(x, "GenomicRatioSet") | is(x, "RGChannelSet")) { 
    if ("XYstats" %in% names(metadata(x))) {
      XY <- as.matrix(metadata(x)$XYstats)
    } else {
      XY <- XYstats(x)
    }
    XY <- XY[, colnames(x)]
    if (is.null(sex) & "sex" %in% names(colData(x))) {
      message("Splitting plot based on `sex`...")
      sex <- x$sex
    } else if (is.null(sex) & "inferred_sex" %in% names(colData(x))) {
      message("Splitting plot based on `inferred_sex`...")
      sex <- x$inferred_sex
    }
  } else { 
    XY <- x
  }
  
  if (!is.null(sex)) {
    Heatmap(XY, name="frac", col=jet, column_split=sex,
            clustering_method_columns='ward.D2', 
            clustering_distance_columns='manhattan', 
            clustering_method_rows='ward.D2', 
            clustering_distance_rows='manhattan',
            ...)
  } else {
    Heatmap(XY, name="frac", col=jet,
            clustering_method_columns='ward.D2', 
            clustering_distance_columns='manhattan', 
            clustering_method_rows='ward.D2', 
            clustering_distance_rows='manhattan',
            ...)

  } 

}
