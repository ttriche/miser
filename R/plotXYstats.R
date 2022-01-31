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
      XY <- as.matrix(metadata(grSet)$XYstats)
    } else {
      XY <- XYstats(grSet)
    }
    XY <- XY[, colnames(x)]
    if (is.null(sex) & "sex" %in% names(colData(x))) {
      message("Splitting plot based on `sex`...")
      sex <- grSet$sex
    } else if (is.null(sex) & "inferred_sex" %in% names(colData(x))) {
      message("Splitting plot based on `inferred_sex`...")
      sex <- grSet$inferred_sex
    }
  } else { 
    XY <- x
  }
  
  if (!is.null(sex)) {
    Heatmap(XY, name="frac", col=jet, column_split=sex, ...)
  } else {
    Heatmap(XY, name="frac", col=jet, ...)
  } 

}
