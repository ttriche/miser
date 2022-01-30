#' dump QC (NAs, flags, XY stats) and (optionally) SNP calls and beta values
#'
#' NA fraction, median CN, control flags, XY stats are saved in stub_QC.csv.gz
#' SNP calls are saved in stub_SNPcalls.csv.gz if snps == TRUE
#' SNP intensities are saved in stub_SNPs.csv.gz if snps == TRUE 
#' Sesamized betas are saved in stub_betas.csv.gz if betas == TRUE 
#'
#' @param grSet       GenomicRatioSet with betas and QC metrics to dump out 
#' @param stub        mandatory prefix for derived files (e.g. stub_QC.csv.gz)
#' @param path        where to save the dumped CSV file[s] (".")
#' @param snps        dump SNPs and SNP calls? (TRUE) 
#' @param betas       dump beta values too? (FALSE) 
#'
#' @return            the status of the final operation
#' 
#' @seealso           SNPcalls
#' 
#' @import            mclust
#' @import            minfi
#' 
#' @export 
dumpQCfiles <- function(grSet, stub=NULL, path=".", snps=TRUE, betas=FALSE) {

  if (is.null(stub) | stub == "") stop("Error: `stub` cannot be empty.")
  if (!is(grSet, "GenomicRatioSet")) stop("Input is not a GenomicRatioSet.")
  grSet <- rename_meta(grSet) # just in case the colnames don't match!

  na_frac <- .NAfrac(grSet)
  XYstats <- .XYstats(grSet)
  medianCN <- .medianCN(grSet)
  flagstats <- metadata(grSet)$control_flagged
  qcfile <- paste0(stub, "_qc.csv")
  qcpath <- file.path(path, qcfile)
  qcmat <- rbind(NA_frac=na_frac, medianCN=medianCN, flagstats, XYstats)
  write.csv(qcmat, qcpath)
  qcgz <- gzip(qcpath)
  message("Wrote QC information to ", qcgz, ".")

  if (snps) {

    SNPcalls <- SNPcalls(grSet)
    SNPcallfile <- paste0(stub, "_SNPcalls.csv")
    SNPcallpath <- file.path(path, SNPcallfile)
    write.csv(SNPcalls, SNPcallpath)
    SNPcallgz <- gzip(SNPcallpath)
    message("Wrote SNP calls to ", SNPcallgz, ".")

    SNPs <- metadata(grSet)$SNPs
    SNPfile <- paste0(stub, "_SNPs.csv")
    SNPpath <- file.path(path, SNPfile)
    write.csv(SNPs, SNPpath)
    SNPgz <- gzip(SNPpath)
    message("Wrote SNP intensities to ", SNPgz, ".")

  }

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
.medianCN <- function(grSet) {
  
  colMedians(getCN(grSet), na.rm=TRUE)

}


# helper fn
.XYstats <- function(grSet) {

  XNA <- is.na(getBeta(subset(grSet, seqnames %in% c("chrX", "X"))))
  rowfrac <- rowSums(XNA) / ncol(XNA)
  masked_X <- sum(rowfrac == 1)
  XNAfrac <- (colSums(XNA) - masked_X) / (nrow(XNA) - masked_X)

  YNA <- is.na(getBeta(subset(grSet, seqnames %in% c("chrY", "Y"))))
  yrowfrac <- rowSums(YNA) / ncol(YNA)
  masked_Y <- sum(yrowfrac == 1)
  YNAfrac <- rep(1, ncol(YNA))
  if (!all(yrowfrac == 1)) {
    YNAfrac <- (colSums(YNA)-masked_Y)/(nrow(YNA)-masked_Y)
  }

  whichX <- rownames(subset(grSet, seqnames %in% c("chrX", "X") & 
                                   rowData(grSet)$IslandStatus != "OpenSea"))
  XBeta <- colMedians(getBeta(grSet[whichX, ]), na.rm=TRUE)
  
  return(rbind(XNAfrac=XNAfrac, YNAfrac=YNAfrac, XBeta=XBeta))

}
