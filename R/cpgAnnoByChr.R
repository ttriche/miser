#' DMRcate doesn't play nicely with HDF5-backed SEs, so I wrote this instead.
#' 
#' TODO: run in parallel, especially for HDF5-backed objects.
#' 
#' @param object      a GenomicRatioSet 
#' @param design      a design matrix 
#' @param coef        which column to fit (default is 2)
#' @param fdr         FDR cutoff (0.05)
#' @param ...         other arguments, to be passed (in principle) to DMRcate
#' 
#' @return        annotated fits for x and its design matrix
#'
#' @import DMRcate
#' @import limma
#' 
#' @export 
cpgAnnoByChr <- function(object, design, coef=2, fdr=0.05, ...) { 

  fit <- function(x) { 
    message("Annotating ", unique(seqnames(x)), "...") 
    res <- list(
      M=eBayes(lmFit(as(getM(x), "matrix"), design)),
      B=eBayes(lmFit(as(getBeta(x), "matrix"), design))
    )
    res$tt <- topTable(res$M, coef=coef, number=nrow(x))
    betatt <- topTable(res$B, coef=coef, number=nrow(x))
    m <- match(rownames(res$tt), rownames(betatt))
    res$tt$betafc <- betatt$logFC[m]
    m <- match(rownames(x), rownames(res$tt))
    res$tt <- res$tt[m,]
    res$anno <- data.frame(ID=rownames(res$tt), stat=res$tt$t, pval=res$tt$P.Value,
                           CHR=seqnames(x), pos=start(x), betafc=res$tt$betafc) 
    return(res$anno)
  }

  # this could easily run in parallel, of course 
  annotated <- do.call(rbind, lapply(split(object, seqnames(object)), fit))
  annotated$indfdr <- p.adjust(annotated$pval, method="BH")
  annotated$is.sig <- (annotated$indfdr < fdr)
  message(sum(annotated$is.sig), " significant ", 
          ifelse(sum(annotated$is.sig) == 1, "locus", "loci"), " found.")
  class(annotated) <- "annot"
  return(annotated)

} 
