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
  if (!"control_flagged" %in% names(metadata(grSet))) stop("No control data!")
  grSet <- rename_meta(grSet) # just in case the colnames don't match!

  qcfile <- paste0(stub, "_qc.csv")
  qcpath <- file.path(path, qcfile)
  qcmat <- rbind(NA_frac=NAfrac(grSet), 
                 medianCN=colMedians(getCN(grSet)),
                 metadata(grSet)$control_flagged,
                 XYstats(grSet))
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
