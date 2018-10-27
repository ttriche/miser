#' sesamize a bunch of IDATs from GEO (or anywhere else with sensible names)
#' 
#' This function processes a directory full of IDATs, somewhat sensibly. 
#'
#' If an RGChannelSet object is provided as the first argument, assume it's 
#' already got a $subject colData column and process as if we'd got to there. 
#'
#' A better-er idea would probably involve running the IDATs individually and
#' then bolting them all together into an HDF5-backed grSet at the last instant.
#' 
#' @param   subjects  the names of each subject, or an rgSet (autodetect)
#' @param   frags     which elements to extract for the array Basename (1:3)
#' @param   annot     look up the titles and characteristics for GSMs? (FALSE)
#' @param   HDF5      EXPERIMENTAL: make the grSet HDF5-backed? (FALSE)
#' @param   ...       more arguments to pass on to sesamize
#'
#' @return            a GenomicRatioSet with metadata(grSet)$SNPs filled out
#' 
#' @import  minfi
#' @import  sesame
#' 
#' @export 
sesamizeGEO <- function(subjects=NULL, frags=1:3, annot=FALSE, HDF5=FALSE, ...){

  if (is(subjects, "RGChannelSet")) {
    stopifnot(all(c("subject","Basename") %in% names(colData(subjects))))
    message("Testing annotations on the first few samples...") 
    tmp <- sesamize(subjects[,1:2]) # will stop() if anything is out of wack
    message("Sesamizing...") 
    res <- sesamask(sesamize(subjects, ...))
  } else { 
    samps <- getSamps(subjects=subjects, frags=frags)
    message("Testing annotations on the first few IDATs...") 
    stopifnot(testIDATs(samps, frags)) 
    message("Reading IDATs into a temporary RedGreenChannelSet...") 
    rgSet <- getRGChannelSet(subjects=subjects, frags=frags)
    message("Sesamizing...") 
    res <- sesamask(sesamize(rgSet, ...))
  }
  if (annot) res <- addTitles(addCharateristics(res))
  if (HDF5) message("HDF5 support isn't set up yet")
  return(res)

}
