#' simple function for plotting control probe failure flaggings
#' 
#' If an object is provided without suitable metadata, this function will fail.
#' See flag_control_failures().
#' 
#' @param x     an object w/!is.null(metadata(x)$control_flagged), or a matrix
#' @param only  only plot samples where there are failures? (FALSE)
#' @param ... additional parameters to pass to ComplexHeatmap::Heatmap
#' 
#' @return a Heatmap
#' 
#' @import ComplexHeatmap
#' @import circlize
#' 
#' @seealso flag_control_failures
#'
#' @export
plotFlaggedControls <- function(x, only=FALSE, ...) {

  if (is(x, "GenomicRatioSet") | is(x, "RGChannelSet")) { 
    if (!"control_flagged" %in% names(metadata(x))) {
      stop("No control_flagged metadate found, cannot proceed.")
    } else {
      flagged <- as.matrix(metadata(x)$control_flagged[, colnames(x)])
    }
  } else { 
    flagged <- as.matrix(x)
  }
  flagged <- flagged[!apply(flagged, 1, function(x) all(is.na(x))), ]
  
  flagcols <- colorRamp2(c(0, 1), c("white", "darkred"))
  all_samples <- ifelse(colSums(flagged, na.rm=TRUE) == 0, 
                        "all metrics pass", "some metrics fail")
  all_probes <- ifelse(rowSums(flagged, na.rm=TRUE) == 0,
                       "all samples pass", "some samples fail")

  if (only) { 
    Heatmap(flagged[, all_samples != "all metrics pass"], 
            name="failed", col=flagcols, 
            clustering_method_columns='ward.D2', 
            clustering_distance_columns='manhattan', 
            row_split=all_probes,
            clustering_method_rows='ward.D2', 
            clustering_distance_rows='manhattan',
            ...) 
  } else { 
    Heatmap(flagged, name="failed", col=flagcols, 
            column_split=all_samples,
            clustering_method_columns='ward.D2', 
            clustering_distance_columns='manhattan', 
            row_split=all_probes,
            clustering_method_rows='ward.D2', 
            clustering_distance_rows='manhattan',
            ...) 
  }

}
