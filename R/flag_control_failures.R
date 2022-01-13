#' get control quality assessments from an object returned by `bactrl`
#' 
#' Apply BeadArray thresholds for quality flagging. Thresholds obtained from 
#' the BeadArray Controls Reporter Software Guide (v00, source 2) and ewastools
#' resource (v1.5, source 3).
#'
#' @param rmat  return matrix from control_metrics(...), N rows by M cols, or
#'              an RGChannelSet (in which case control_metrics(x) will be run)
#' @param dft   optional threshold data.frame; defaults will be used if null 
#' 
#' @return      an NxM matrix of threshold assessments, 0 = pass, 1 = fail
#' 
#' @export 
flag_control_failures <- function(rmat, dft=NULL) {

  # if it's an RGChannelSet, generate metrics for it 
  if (is(rmat, "RGChannelSet")) rmat <- control_metrics(rmat, dft=dft)

  fails <- rmat
  if (is.null(dft)) dft <- .get_dft()
  for (i in colnames(rmat)[-1]) fails[, i] <- as.numeric(rmat[, i] < dft[, i])
  return(fails)

}
