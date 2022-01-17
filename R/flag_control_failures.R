#' Get control QA from an object returned by `bactrl` or from an RGChannelSet
#' 
#' Apply BeadArray thresholds for quality flagging. Thresholds obtained from 
#' the BeadArray Controls Reporter Software Guide (v00, source 2) and ewastools
#' resource (v1.5, source 3). If an RGChannelSet is passed in, control_metrics()
#' will be run on the rgSet and its detected annotation, then failures flagged.
#'
#' @param rmat      return matrix from control_metrics(), or RGChannelSet 
#' @param dft       optional threshold data.frame; defaults will be used if null
#' @param platform  optional platform indicator; will detect if rgSet ("EPIC") 
#' 
#' @return      an NxM matrix of threshold assessments, 0 = pass, 1 = fail
#'
#' @examples
#' 
#' if (exists("MPAL_rgSet")) {
#' 
#'   # TARGET MPAL data, from Alexander et al, Nature 2018
#'   flagged <- flag_control_failures(MPAL_rgSet)
#'   all_samples_passed <- ifelse(colSums(flagged) == 0, 1, 0)
#'   all_probes_passed <- ifelse(rowSums(flagged) == 0, 1, 0)
#' 
#'   library(circlize)
#'   flagcols <- colorRamp2(c(0, 1), c("white", "darkred"))
#' 
#'   library(ComplexHeatmap)
#'   Heatmap(t(flagged), name="failed", col=flagcols, column_names_side="top",
#'           column_split=all_probes_passed, column_names_gp=gpar(fontsize=6), 
#'           row_split=all_samples_passed, row_names_side="left")
#'
#' } 
#'   
#' @details
#' Every single EPIC IDAT that we have seen fails nonpolymorphic.grn and 99.5% 
#' of EPIC arrays fail bisulfite.conv.I.grn, so we mask these (for now) on EPIC.
#' Zhou, Laird, and Shen (NAR 2017) document additional reasons why bisulfite
#' conversion control probes on this platform are probably not informative, not
#' least because poor bisulfite conversion will typically result in poor
#' hybridization due to non-target cytosine bases in most probe oligos. In 
#' conclusion, these two control probe sets seem to be uninformative on EPIC,
#' so they will be masked as NA if an EPIC rgSet is passed in. 
#'
#' @seealso         control_metrics
#'
#' @export 
flag_control_failures <- function(rmat, dft=NULL, platform="epic") {

  # if it's an RGChannelSet, generate metrics for it 
  if (is(rmat, "RGChannelSet")) {
    platform <- tolower(sub("IlluminaHumanMethylation", "", 
                            annotation(rmat)["array"]))
    rmat <- control_metrics(rmat, dft=dft)
  }

  fails <- rmat
  if (is.null(dft)) dft <- .get_dft()
  for (i in colnames(rmat)) fails[, i] <- as.numeric(rmat[, i] < dft[, i])
  if (platform == "epic") {
    message("Masking dubious results on EPIC array control probes...")
    is.na(fails[, c("bisulfite.conv.I.grn","nonpolymorphic.grn")]) <- TRUE
  }

  return(fails)

}
