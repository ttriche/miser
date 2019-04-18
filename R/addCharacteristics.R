#' get the characteristics of a bunch of GSM entries
#' 
#' like it says on the tin. this is usually for filling out pData/colData
#' 
#' @param   x         either a grSet from sesamizeGEO, or a bunch of GSMs
#' @param   column    name of the column holding the GSMs ("subject") 
#' @param   titles    tack on the titles as well? (TRUE)
#' @param   cachePath where to cache the GEOmetadb sqlite file (tempdir())
#'
#' @return            a grSet, perhaps with additional colData, or a data.frame
#' 
#' @import  GEOmetadb
#' @import  RSQLite
#' 
#' @export 
addCharacteristics <- function(x, column="subject",titles=TRUE,cachePath=NULL) {
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
  if (!all(GSMs %in% res$gsm)) stop("Missing records! Aborting.")
  toParse <- res[GSMs, "characteristics_ch1"]
  parsed <- do.call(rbind, lapply(strsplit(toParse,";\t"), elts,sep=": ",elt=2))
  colnames(parsed) <- sapply(strsplit(toParse[1], ";\t")[[1]], elts, sep=":")
  rownames(parsed) <- GSMs
  parsed <- as.data.frame(parsed)
  if (titles) parsed$title <- addTitles(GSMs)
  parsed$gsm <- GSMs
  if (is(x, "GenomicRatioSet")) {
    for (i in names(parsed)) colData(x)[, i] <- parsed[,i]
    return(x)
  } else {
    return(parsed)
  }
}
