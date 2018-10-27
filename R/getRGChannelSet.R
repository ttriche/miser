#' assemble a suitable RGChannelSet for sesamizing
#' 
#' Note: this should probably go away in favor of signalSet parallel processing
#' 
#' @param   subjects  the names of each subject, or NULL to autodetect (NULL)
#' @param   frags     which elements to extract for the array Basename (1:3)
#' 
#' @return            an RGChannelSet with pData $subject and $Basename filled 
#' 
#' @import  minfi
#' @import  sesame
#' 
#' @export 
getRGChannelSet <- function(subjects=NULL, frags=1:3) { 
  samps <- getSamps(subjects=subjects, frags=frags)
  rgSet <- read.metharray.exp(base=".", targets=samps, verbose=TRUE)
  sampleNames(rgSet) <- rgSet$subject
  return(rgSet) 
} 
