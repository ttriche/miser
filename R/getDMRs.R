#' use DMRcate to get DMRs (with less hassle) 
#' 
#' @param x               a GenomicRatioSet 
#' @param design          a design matrix (NULL)
#' @param dropXY          drop sex chromosomes? (TRUE, usually a good idea)
#' @param impute          impute NAs if > 0.5 of samples are non-NA? (TRUE)
#' @param contrasts       apply contrasts? (FALSE)
#' @param cont.matrix     a contrast matrix (NULL) 
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
getDMRs <- function(x, design=NULL, dropXY=TRUE, impute=TRUE, contrasts=FALSE, 
                    cont.matrix=NULL, coef=2, fdr=.05, betacutoff=.1, DMLs=TRUE,
                    ...) {

  if (dropXY) { 
    message("Dropping sex chromosomes (dropXY is set to TRUE)...") 
    x <- keepSeqlevels(x, paste0("chr", 1:22), pruning.mode="coarse")
  }
  if (any(is.na(getBeta(x))) & impute) {
    message("Imputing NAs...") 
    x <- fixNAs(x)
  }
  message("Annotating individual CpGs...")
  annot <- cpg.annotate(datatype="array", x, coef=coef, 
                        design=design, contrasts=contrasts, 
                        cont.matrix=cont.matrix, fdr=fdr)
  message("Demarcating significant regions...")
  res <- dmrcate(annot, betacutoff=betacutoff, ...) 
  res$granges <- extractRanges(res, genome=unique(genome(x)))
  if (DMLs) res$DMLs <- annot
  return(res)

}
