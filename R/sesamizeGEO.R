#' sesamize a bunch of IDATs from GEO (or anywhere else with sensible names)
#' 
#' I got tired of typing the same crap to process IDATs.
#' This function processes a directory full of IDATs, but sensibly. 
#' 
#' @param   subjects  the names of each subject, if known (default: autodetect)
#' @param   elts      which elements to extract for the array Basename (1:3)
#' @param   mask      add rowData(grSet)$mask using sesameData? (TRUE) 
#' @param   ...       more arguments to pass on to sesamize
#'
#' @return            a GenomicRatioSet with metadata(grSet)$SNPs filled out
#' 
#' @import  minfi
#' @import  sesame
#' 
#' @export 
sesamizeGEO <- function(subjects=NULL, elts=c(1,2,3), mask=TRUE, ...) { 
  samps <- data.frame(Basename=unique(elts(list.files(patt="*idat*"), z=elts)))
  if (is.null(subjects)) {
    samps$subject <- elts(samps$Basename)
  } else { 
    stopifnot(identical(names(subjects), samps$Basename))
  }
  message("Reading IDATs into temporary RedGreenChannelSet...") 
  rgSet <- read.metharray.exp(base=".", targets=samps, verbose=TRUE)
  sampleNames(rgSet) <- rgSet$subject
  message("Sesamizing...") 
  res <- sesamize(rgSet, ...)
  if (mask) {
    message("Adding rowData(grSet)$mask...")
    return(sesamask(res))
  } else { 
    return(res)
  } 
}
