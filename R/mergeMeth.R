#' helper function to merge grSets
#' 
#' @param x       a GenomicRatioSet
#' @param y       another GenomicRatioSet
#' @param verbose be verbose? (TRUE) 
#' 
#' @return        a GenomicRatioSet combining the two
#' 
#' @import        S4Vectors
#' @import        minfi
#' 
#' @examples
#' # merged <- Reduce(mergeMeth, list(grs1, grs2, grs3))
#' 
#' @export 
mergeMeth <- function(x, y, verbose=TRUE) { 

  stopifnot(is(x, "GenomicRatioSet"))
  stopifnot(is(y, "GenomicRatioSet"))

  # assays
  matchingAssays <- intersect(assayNames(x), assayNames(y))
  if (length(matchingAssays) < 1) stop("No matching assays found.")
  if (verbose) message(length(matchingAssays), " assays matched.")
  assays(x) <- assays(x)[matchingAssays]
  assays(y) <- assays(x)[matchingAssays]

  # probes 
  matchingRows <- intersect(rownames(x), rownames(y))
  if (length(matchingRows) < 1) stop("No matching probes found.")
  if (verbose) message(length(matchingRows), " probes matched.")
  # x <- x[matchingRows,] # don't bother
  # y <- y[matchingRows,] # don't bother 

  # covariates
  matchingColData <- intersect(names(colData(x)), names(colData(y)))
  if (length(matchingColData) < 1) warning("No matching colData columns found.")
  if (verbose) message(length(matchingColData), " colData columns matched.")
  colData(x) <- colData(x)[, matchingColData]
  colData(y) <- colData(y)[, matchingColData]

  # if no SNPs in metadata, this is the last step
  res <- cbind(x[matchingRows,], y[matchingRows,])
 
  # metadata/SNP weirdness fix 
  if ("SNPs" %in% names(metadata(x)) & "SNPs" %in% names(metadata(y))) {
    SNPs <- cbind(metadata(x)$SNPs[, colnames(x)], 
                  metadata(y)$SNPs[, colnames(y)])
    metadata(res) <- metadata(res)[names(metadata(res)) != "SNPs"]
    metadata(res)$SNPs <- SNPs # yes this is bizarre, yes it is needed
  }

  return(res)

}
