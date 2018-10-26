#' use DMRcate to get DMRs (with less hassle) 
#' 
#' @param x               a GenomicRatioSet 
#' @param design          a design matrix (NULL)
#' @param contrasts       apply contrasts? (FALSE)
#' @param cont.matrix     a contrast matrix (NULL) 
#' @param coef            which column to fit (or "all") (default is 2)
#' @param fdr             FDR cutoff (0.05)
#' @param betacutoff      DMRs must have at least this maxbetaFC (0.1)
#' @param keep.raw.DMLs   keep differentially methylated loci for plots? (TRUE)
#' @param ...             other arguments to pass to dmrcate
#' 
#' @return                a dmrcate.output object
#' 
#' @import DMRcate
#'
#' @export
getDMRs <- function(x, design=NULL, contrasts=FALSE, cont.matrix=NULL, coef=2, 
                    fdr=.05, betacutoff=.1, keep.raw.DMLs=TRUE, ...) {

  message("Annotating individual CpGs...")
  annot <- cpg.annotate(datatype="array", x, coef=coef, 
                        design=design, contrasts=contrasts, 
                        cont.matrix=cont.matrix, fdr=fdr)
  message("Demarcating significant regions...")
  res <- dmrcate(annot, betacutoff=betacutoff, ...) 
  res$granges <- extractRanges(res, genome=unique(genome(x)))
  if (keep.raw.DMLs) res$DMLs <- annot
  return(res)

}
