#' convenience function for merging large RGChannelSets
#' 
#' Merges metadata, drops from each RGChannelSet, merges RGChannelSets, re-adds
#' RGChannelSets' metadata must have identical names, and platform must match.
#' 
#' @param rgSet1  the first RGChannelSet
#' @param rgSet2  the second RGChannelSet
#' 
#' @return        the result of cbind(rgSet1, rgSet2) with sensible metadata 
#' 
#' @import        minfi
#' 
#' @export
merge_rgSets <- function(rgSet1, rgSet2) { 

  stopifnot(is(rgSet1, "RGChannelSet"))
  stopifnot(is(rgSet2, "RGChannelSet"))
  stopifnot(identical(annotation(rgSet1), annotation(rgSet2)))
  shared <- intersect(rownames(rgSet1), rownames(rgSet2))
  if (length(shared) == 0) {
    stop("No shared rows between rgSet1 and rgSet2, cannot proceed.")
  } else if (!identical(rownames(rgSet1), rownames(rgSet2))) {
    warning("Rownames differ. Keeping ", length(shared), " shared rows.")
  }
  stopifnot(identical(names(colData(rgSet1)), names(colData(rgSet2))))
  if (any(duplicated(c(colnames(rgSet1), colnames(rgSet2))))) {
    stop("Duplicate column names in merged RGChannelSet, cannot proceed.")
  }

  md <- merge_metadata(rgSet1, rgSet2)
  metadata(rgSet1) <- list()
  metadata(rgSet2) <- list() 

  message("Merging ", ncol(rgSet1), " samples and ", ncol(rgSet2), " samples.")
  rgSet <- cbind(rgSet1[shared, ], rgSet2[shared, ])
  metadata(rgSet) <- md
  message("Done. Merged RGChannelSet contains ", ncol(rgSet), " samples.")
  return(rgSet)

}
