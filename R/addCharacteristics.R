#' get the characteristics of a bunch of GSM entries
#' 
#' like it says on the tin. this is usually for filling out pData/colData
#' 
#' @param   x         an RGChannelSet or GenomicRatioSet with $subject, or GSMs
#' @param   column    name of the column holding the GSMs ("subject") 
#' @param   titles    tack on the titles as well? (TRUE)
#' 
#' @details           this function rolls through the GSMs and grabs @header
#'
#' @return            an rg/grSet, perhaps with more colData, or a data.frame
#' 
#' @import  GEOquery
#' 
#' @export 
addCharacteristics <- function(x, column="subject",titles=TRUE,cachePath=NULL) {

  if (is(x, "GenomicRatioSet") | is(x, "RGChannelSet")) {
    if (!column %in% names(colData(x))) {
      stop("You need a column named ", column, " in your colData to run this")
    } else { 
      GSMs <- x$subject
    }
  } else { 
    GSMs <- x
  }


  if (!all(GSMs %in% res$gsm)) stop("Missing records! Aborting.")
  
  toParse <- res[GSMs, "characteristics_ch1"]
  parsed <- do.call(rbind, lapply(strsplit(toParse,";\t"), elts,sep=": ",elt=2))
  colnames(parsed) <- sapply(strsplit(toParse[1], ";\t")[[1]], elts, sep=":")
  rownames(parsed) <- GSMs
  parsed <- as.data.frame(parsed)
  if (titles) parsed$title <- addTitles(GSMs)
  parsed$gsm <- GSMs

  if (is(x, "GenomicRatioSet") | is(x, "RGChannelSet")) {
    for (i in names(parsed)) colData(x)[, i] <- parsed[,i]
    return(x)
  } else {
    return(parsed)
  }

}
