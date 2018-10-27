#' helper function to read in samples and Basename columns
#' 
#' @param   subjects  the names of each subject, or NULL to autodetect (NULL)
#' @param   frags     which elements to extract for the array Basename (1:3)
#' 
#' @return            a data.frame
#' 
#' @export 
getSamps <- function(subjects=NULL, frags=1:3) { 
  basenames <- unique(elts(list.files(patt="*idat*"), z=frags))
  samps <- data.frame(Basename=basenames)
  if (is.null(subjects)) {
    samps$subject <- elts(samps$Basename)
  } else { 
    stopifnot(identical(names(subjects), samps$Basename))
  }
  return(samps)
}
