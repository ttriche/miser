#' assemble a suitable RGChannelSet for sesamizing
#' 
#' Note: this should probably go away in favor of signalSet parallel processing
#' 
#' @param   subjects  the names of each subject, or NULL to autodetect (NULL)
#' @param   frags     which elements to extract for the array Basename (1:3)
#' @param   samps     targets data.frame (overrides subjects and frags) (NULL)
#' @param   ...               additional arguments to pass to read.metharray()
#' 
#' @return            an RGChannelSet with pData $subject and $Basename filled 
#' 
#' @import  minfi
#' @import  sesame
#' 
#' @export 
getRGChannelSet <- function(subjects=NULL, frags=1:3, samps=NULL, ...) { 
  if (is.null(samps)) samps <- getSamps(subjects=subjects, frags=frags)
  stopifnot(all(c("Basename","subject") %in% names(samps)))
  rgSet <- read.metharray.exp(base=".", targets=samps, verbose=TRUE, ...)
  sampleNames(rgSet) <- rgSet$subject
  return(rgSet) 
} 
