#' convenience function to rename columns of *Set metadata matrices
#' 
#' @param   x   the SummarizedExperiment-derived thing with metadata to rename
#' 
#' @return      the same, but with metadata columns renamed to colnames(x)
#' 
#' @import      SummarizedExperiment
#'
#' @export
rename_meta <- function(x) {
  mncols <- vapply(metadata(x), ncol, integer(1))
  if (length(mncols) == 0) {
    message("No metadata found. Returning unchanged.")
  } else if (!all(mncols == ncol(x))) {
    for (i in names(mncols)[which(mncols != ncol(x))]) {
      message("ncol(metadata(x)[[", i, "]]) != ncol(x).")
    }
    stop("Fix your metadata before running this function!")
  } else { 
    for (i in names(mncols)) {
      colnames(metadata(x)[[i]]) <- colnames(x)
    }
  }
  return(x) 
}
