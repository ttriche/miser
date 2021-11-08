#' convenience function; checks to see if an IDAT file was downloaded already
#' 
#' @param IDAT  the path to an IDAT (e.g. getGEO(GSE)[[1]]$supplemental_file[1])
#' @param path  where to look for it or store it locally (".")
#' 
#' @return      the local file path, if present or if downloaded successfully
#' 
#' @import      utils
#'
#' @export
getIDAT <- function(IDAT, path=".") {

  if (length(IDAT) > 1) return(vapply(IDAT, getIDAT, character(1)))

  fname <- basename(IDAT)
  if (fname %in% list.files(path)) {
    message("Found ", fname, " locally; skipping download.")
  } else {
    message("Downloading ", fname, " from ", IDAT, "...", appendLF=FALSE)
    download.file(IDAT, file.path(path, fname))
    message("done.")
  }
  
  if (path == ".") { 
    return(fname)
  } else { 
    return(file.path(path, fname))
  }

}
