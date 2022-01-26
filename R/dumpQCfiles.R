#' dump QC (NAs, flags) and SNP/XY covariates CSVs (optionally, betas as well)
#'
#' SNPs and sex are saved in stub_SNPXY.csv.gz
#' NA fraction and flagged controls are saved in stub_qc.csv.gz
#' sesamized betas (including NA-masked values) are saved in stub_betas.csv.gz
#' (but only if betas == TRUE! this tends to be an enormous file, obviously) 
#'
#' @param grSet       GenomicRatioSet with betas and QC metrics to dump out 
#' @param stub        mandatory file prefix (_qc, _SNPXY, _betas.csv.gz)
#' @param path        where to save the dumped CSV file[s] (".")
#' @param betas       dump beta values too? (FALSE) 
#'
#' @return a GenomicRatioSet (or an rgSet in case of failure)
#' 
#' @import minfi
#' 
#' @export 
dumpQCfiles <- function(grSet, stub, path=".", betas=FALSE) {

  if (!is(grSet, "GenomicRatioSet")) stop("Input is not a GenomicRatioSet.")
  grSet <- rename_meta(grSet) # just in case the colnames don't match!

  na_frac <- .NAfrac(grSet)
  qcfile <- paste0(stub, "_qc.csv")
  qcpath <- file.path(path, qcfile)
  qcmat <- rbind(NA_frac=na_frac, metadata(grSet)$control_flagged)
  write.csv(qcmat, qcpath)
  qcgz <- gzip(qcpath)
  message("Wrote QC information to ", qcgz, ".")

  snps_xy <- .SNPXY(grSet)
  snpxyfile <- paste0(stub, "_SNPXY.csv")
  snpxypath <- file.path(path, snpxyfile)
  write.csv(snps_xy, snpxypath)
  snpxygz <- gzip(snpxypath)
  message("Wrote SNPXY information to ", snpxygz, ".")

  if (betas) { 
    
    betafile <- paste0(stub, "_betas.csv")
    betapath <- file.path(path, betafile)
    write.csv(getBeta(grSet), betapath)
    betagz <- gzip(betapath)
    message("Wrote betas to ", betagz, ".")

  }

}


# helper fn
.NAfrac <- function(grSet) {

  na_mat <- is.na(getBeta(grSet))
  rowfrac <- rowSums(na_mat) / ncol(na_mat)
  masked_loci <- sum(rowfrac == 1)
  message(masked_loci, " are masked in all samples; omitting from NA tallies.")
  na_frac <- (colSums(na_mat) - masked_loci) / (nrow(na_mat) - masked_loci)
  worstsample <- colnames(na_mat)[which.max(na_frac)]
  message("Failed probes: median of ", round(median(na_frac * 100), 1), "%, ", 
          "maximum of ", round(max(na_frac * 100), 1), "% ",
          "(in ", worstsample, ")")
  return(na_frac) 

}


# helper fn
.SNPXY <- function(grSet) {

  SNPs <- metadata(grSet)$SNPs 

  XNA <- is.na(getBeta(subset(grSet, seqnames %in% c("chrX", "X"))))
  rowfrac <- rowSums(XNA) / ncol(XNA)
  masked_X <- sum(rowfrac == 1)
  XNAfrac <- (colSums(XNA) - masked_X) / (nrow(XNA) - masked_X)


  YNA <- is.na(getBeta(subset(grSet, seqnames %in% c("chrY", "Y"))))
  yrowfrac <- rowSums(YNA) / ncol(YNA)
  masked_Y <- sum(yrowfrac == 1)
  if (all(yrowfrac == 1)) {
    YNAfrac <- rep(1, ncol(YNA))
  } else { 
    YNAfrac <- (colSums(YNA) - masked_Y) / (nrow(YNA) - masked_Y)
  }

  whichX <- rownames(subset(grSet, 
                            seqnames %in% c("chrX", "X") & 
                              rowData(grSet)$IslandStatus != "OpenSea"))
  XBeta <- colMeans(getBeta(grSet[whichX, ]), na.rm=TRUE)
  rbind(SNPs, XNAfrac=XNAfrac, YNAfrac=YNAfrac, XBeta=XBeta)

}
