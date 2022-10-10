#' convenience function for merging large SummarizedExperiment-like objects
#' 
#' The two objects' metadata must have identical names,
#' and each metadata items' rownames must be identical.
#' 
#' @param x       the first object 
#' @param y       the second object 
#' 
#' @return        a list() with the merged metadata
#' 
#' @export
merge_metadata <- function(x, y) {

  rename_meta(x)
  md1 <- metadata(x)
  stopifnot(identical(colnames(md1[[1]]), colnames(x)))

  rename_meta(y)
  md2 <- metadata(y) 
  stopifnot(identical(colnames(md2[[1]]), colnames(y)))

  stopifnot(identical(names(md1), names(md2)))
  md <- list() 
  for (i in names(md1)) {
    message("Merging metadata item `", i, "`...", appendLF=FALSE)
    stopifnot(identical(rownames(md1[[i]]), rownames(md2[[i]])))
    md[[i]] <- cbind(md1[[i]], md2[[i]])
    message("Done.")
  }
  return(md)

}
