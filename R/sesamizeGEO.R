#' sesamize a bunch of IDATs from GEO (or anywhere else with sensible names)
#' 
#' This function processes a directory full of IDATs, but sensibly. 
#' 
#' @param   subjects  the names of each subject, if known (default: autodetect)
#' @param   elts      which elements to extract for the array Basename (1:3)
#' @param   mask      add rowData(grSet)$mask using sesameData? (TRUE) 
#' @param   titles    look up the titles for GSM entries from GEO? (FALSE) 
#' @param   ...       more arguments to pass on to sesamize
#'
#' @return            a GenomicRatioSet with metadata(grSet)$SNPs filled out
#' 
#' @import  minfi
#' @import  sesame
#' 
#' @export 
sesamizeGEO <- function(subjects=NULL, elts=1:3, mask=TRUE, titles=FALSE, ...) {
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
  if (titles) res <- titleGEO(res)
  if (mask) res <- sesamask(res)
  return(res)
}
