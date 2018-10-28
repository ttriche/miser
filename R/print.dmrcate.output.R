#' since I use DMRcate a lot, I wrote this to avoid borking the command line
#' 
#' @param object  an object of class dmrcate.output
#' 
#' @return        nothing, just cat the object's various properties
#' 
#' @import        DMRcate
#'
#' @export
print.dmrcate.output <- function(object) {
  objName <- deparse(match.call()$object)
  cat("Object of class", class(object), "\n\n")
  if ("DMLs" %in% names(object)) cat("  Raw DMLs saved in object$DMLs.\n")
  if ("DMRs" %in% names(object)) cat("  GRanges of DMRs in object$DMRs.\n")
  cat("  Input data:", nrow(object$input), "features\n")
  cat("  Stouffer p:", paste(format(range(object$results$Stouffer), 
                                    sci=TRUE, digits=2), collapse=" to "), "\n")
  cat("     Results:", nrow(object$results), "DMRs\n")
}
