#' simple function for plotting control probe failure flaggings
#' 
#' The data is drawn from metadata(grSet)$control_flagged; if not present, 
#' the function will fail with an error message. See flag_control_failures().
#' 
#' @param grSet GenomicRatioSet where !is.null(metadata(grSet)$control_flagged)
#' @param ...   additional parameters to pass to ComplexHeatmap::Heatmap
#' 
#' @return a Heatmap
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @export
plotFlaggedControls <- function(grSet, ...) {

  stopifnot("control_flagged" %in% names(metadata(grSet)))
  flagged <- metadata(grSet)$control_flagged
  flagged <- flagged[!apply(flagged, 1, function(x) any(is.na(x))), ]
  flagcols <- colorRamp2(c(0, 1), c("white", "darkred"))
  all_samples <- ifelse(colSums(flagged, na.rm=TRUE) == 0, 
                        "all metrics pass", "some metrics fail")
  all_probes <- ifelse(rowSums(flagged, na.rm=TRUE) == 0,
                       "all samples pass", "some samples fail")

  Heatmap(flagged, name="failed", col=flagcols, column_names_side="top",
          column_split=all_samples, row_split=all_probes, row_names_side="left")

}
