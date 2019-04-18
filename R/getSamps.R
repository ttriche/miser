#' helper function to read in samples and Basename columns
#' 
#' @param   subjects  the names of each subject, or NULL to autodetect (NULL)
#' @param   frags     which elements to extract for the array Basename (1:3)
#' 
#' @return            a data.frame
#' 
#' @export 
getSamps <- function(subjects=NULL, frags=1:3) { 
  # short circuit all of this if the subjects are a list of IDATs:
  if (all(grepl("idat", ignore.case=TRUE, subjects))) { 
    # a list of IDATs; don't screw around with anything else
    samps <- data.frame(Basename=unique(elts(subjects, elt=frags)),
                        subject=unique(elts(subjects)))
  } else { 
    basenames <- unique(elts(list.files(patt="*idat*"), elt=frags))
    samps <- data.frame(Basename=basenames)
    if (is.null(subjects)) {
      samps$subject <- elts(samps$Basename)
    } else { 
      stopifnot(identical(names(subjects), samps$Basename))
    }
  }
  return(samps)
}
