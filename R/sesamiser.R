#' like getGEO, if getGEO worked directly on IDATs and masked off deletions 
#'
#' Note 0: this function will download the data with getGEO first (!)
#' Note 1: you may fare better with recountmethylation, check there first.
#' Note 2: this would work just as well off the moby recountmethylation rgset.
#' Note 3: this function will exit with an error if no IDATs are available. 
#' Note 4: if there is one supplementary file column, we will need to split it.
#' Note 5: there are some mysterious new errors from sesame::sesamize() lately.
#'
#' @param GSE               which GSE to get
#' @param element           which element of the getGEO result to work on? (1)
#' @param suppcols          cols ("supplementary_file", "supplementary_file.1")
#' @param ...               additional arguments to pass to read.metharray()
#' 
#' @return                  a sesamized GenomicRatioSet
#'
#' @import GEOquery
#' @import sesame
#' @import minfi
#'
#' @export
sesamiser <- function(GSE, element=1, suppcols=c("supplementary_file","supplementary_file.1"), path=".", ...){

  tmp <- as(getGEO(GSE)[[element]], "SummarizedExperiment") 
  if (! all(suppcols %in% names(colData(tmp)))) {
    message("Columns ", paste(suppcols, collapse=", "), " are required.")
    for (suppcol in setdiff(suppcols, names(colData(tmp)))) {
      message("This GSE does not have a column `", suppcol, "`.")
    } 
    stop("Exiting!")
  } else {
    covs <- as(colData(tmp), "data.frame")
    for (suppcol in suppcols) IDATs <- getIDAT(covs[, suppcol])
    names(IDATs) <- covs$geo_accession
    covs <- cbind(covs, getSamps(IDATs))
    message("Reading in IDATs...")
    rgSet <- getRGChannelSet(samps=covs, ...)
    message("Sesamizing...")
    res <- sesame::sesamize(rgSet)
    message("Adding probe mask...")
    res <- sesamask(res)
    message("Done.")
    return(res)
  }

}
