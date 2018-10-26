#' get the characteristics of a bunch of GSM entries
#' 
#' like it says on the tin. this is usually for filling out pData/colData
#' 
#' @param   x         either a grSet from sesamizeGEO, or a bunch of GSMs
#' @param   column    name of the column holding the GSMs ("subject") 
#' @param   titles    tack on the titles as well? (FALSE)
#' @param   cachePath where to cache the GEOmetadb sqlite file (tempdir())
#'
#' @return            whatever is in characteristics_ch1, parsed 
#' 
#' @import  GEOmetadb
#' @import  RSQLite
#' 
#' @export 
characteristics <- function(x, column="subject", titles=FALSE, cachePath=NULL) {
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
  query <- paste0("SELECT gsm, characteristics_ch1 FROM gsm WHERE gsm IN ('",
                  paste(GSMs, collapse="','"), "')")
  res <- dbGetQuery(con, query)
  rownames(res) <- res$gsm
  toParse <- res[GSMs, "characteristics_ch1"]
  parsed <- do.call(rbind, lapply(strsplit(toParse, ";\t"), elts, y=": ", 2))
  colnames(parsed) <- sapply(strsplit(toParse[1], ";\t")[[1]], elts, y=":")
  rownames(parsed) <- GSMs
  parsed <- as.data.frame(parsed)
  if (titles) parsed$title <- titles(GSMs)
  return(parsed)
}
