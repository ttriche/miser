#' get the titles of a bunch of GSM entries
#' 
#' like it says on the tin. this is usually for labeling purposes. 
#' 
#' @param   x         either a grSet from sesamizeGEO, or a bunch of GSMs
#' @param   column    name of the column holding the GSMs ("subject") 
#' @param   cachePath where to cache the GEOmetadb sqlite file (tempdir())
#'
#' @return            a retitled grSet (if x is a grSet) or a bunch of titles
#' 
#' @import  GEOmetadb
#' @import  RSQLite
#' 
#' @export 
titles <- function(x, column="subject", cachePath=NULL) { 
  if (is(x, "GenomicRatioSet")) {
    if (!column %in% names(colData(x))) {
      stop("You need a column named ", column, " in your colData to run this")
    }
    GSMs <- x$subject
  } else { 
    GSMs <- x
  }
  if (is.null(cachePath)) cachePath <- tempdir()
  cacheFile <- paste(cachePath, "GEOmetadb.sqlite", sep="/")
  if (!file.exists(cacheFile)) {
    message("Caching GEOmetadb database...") 
    getSQLiteFile(destdir=cachePath)
  }
  con <- dbConnect(RSQLite::SQLite(), cacheFile)
  query <- paste0("SELECT gsm, title FROM gsm WHERE gsm IN ('",
                  paste(GSMs, collapse="','"), "')")
  res <- dbGetQuery(con, query)
  rownames(res) <- res$gsm
  res <- res[GSMs,] # enforce
  if (is(x, "GenomicRatioSet")) {
    colnames(x) <- x$title <- res$title
    return(x)
  } else { 
    return(res$title)
  }
}
