#' use DMRcate to get DMRs (with less hassle) 
#' 
#' TODO: run in parallel, especially for HDF5 
#' 
#' @param x               a GenomicRatioSet 
#' @param design          a design matrix (NULL)
#' @param dropXY          drop sex chromosomes? (TRUE, usually a good idea)
#' @param impute          impute NAs (TRUE, but beware, this WILL hang on HDF5)
#' @param coef            which column to fit (or "all") (default is 2)
#' @param fdr             FDR cutoff (0.05)
#' @param betacutoff      DMRs must have at least this maxbetaFC (0.1)
#' @param DMLs            keep raw differentially methylated loci to plot (TRUE)
#' @param ...             other arguments to pass to dmrcate
#' 
#' @return                a dmrcate.output object
#' 
#' @import DMRcate
#'
#' @export
getDMRs <- function(x, design=NULL, dropXY=TRUE, impute=TRUE, coef=2, fdr=.05, 
                    betacutoff=.1, DMLs=TRUE, ...) {

  message("Note: this is a simplified, differential-only version of DMRcate.")
  message("DMRcate is capable of much more involved analyses. Read its manual.")

  if (dropXY) {
    message("Dropping sex chromosomes if present (dropXY is set to TRUE)...") 
    x <- keepSeqlevels(x, paste0("chr", 1:22), pruning.mode="coarse")
  }
  
  # try and avoid hanging with HDF5-backed SummarizedExperiments
  if (impute) x <- fixNAs(x) # the alternative is to fail...

  message("Annotating individual CpGs...")
  annot <- cpgAnnoByChr(subset(object, rowData(object)$mask == FALSE), 
                        design=design, coef=coef, fdr=fdr)
  message("Demarcating significant regions...")
  res <- dmrcate(annot, betacutoff=betacutoff, ...) 
  if (DMLs) res$DMLs <- annot
  return(res)

}
