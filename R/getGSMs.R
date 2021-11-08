#' fetch supplemental data (IDATs) for a given dataset 
#' 
#' There's probably a better way to do this, but we haven't found it. 
#' IDAT files will be dumped into the working directory, so setwd() first.
#'
#' @param GSMs        a vector of GEO sample IDs 
#' @param cachePath   where the cached GEOmetadb SQLite file lives (tempdir())
#' 
#' @return            a list of supplementary filenames, usually 2 per GSM
#' 
#' @import GEOmetadb
#' @import curl 
#'
#' @export
#'
getGSMs <- function(GSMs, cachePath=NULL) {

  # boilerplate, refactor this 
  if (is.null(cachePath)) cachePath <- tempdir()
  cacheFile <- paste(cachePath, "GEOmetadb.sqlite", sep="/")
  if (!file.exists(cacheFile)) {
    message("Caching GEOmetadb database...") 
    getSQLiteFile(destdir=cachePath)
  }
  con <- dbConnect(RSQLite::SQLite(), cacheFile)
  query <- paste0("SELECT gsm, supplementary_file FROM gsm WHERE gsm IN ('",
                  paste(GSMs, collapse="','"), "')")
  res <- dbGetQuery(con, query)
  IDATs <- unlist(strsplit(res$supplementary_file, ";\t"))
  numfiles <- length(IDATs)
  numGSMs <- length(GSMs)
  if ((numfiles/numGSMs) < 2) message("Some samples have no IDATs to fetch.") 
  fnames <- sapply(IDATs, getIDAT)
  invisible(fnames)

}
