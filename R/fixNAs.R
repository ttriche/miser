#' idiotically simple masking and imputation for sesamized methylation data 
#'
#' TODO: parallelize for HDF5-backed instances instead of running blockwise 
#'
#' @param x       a GenomicRatioSet
#' @param na      probes with this proportion NA (or higher) are dropped (0.5)
#' @param sqz     a squeeze factor to regularize imputation (1-1e-6 by default) 
#' @param verbose be verbose (eg. for debugging HDF5-backed imputation)? (TRUE)
#' 
#' @return    a GenomicRatioSet with masked probes dropped and NAs imputed 
#' 
#' @import    minfi 
#' @import    impute
#' @import    DelayedMatrixStats
#' 
#' @export
fixNAs <- function(x, na=0.5, sqz=(1 - 1e-6), verbose=TRUE) { 

  stopifnot(is(x, "GenomicRatioSet")) 
  stopifnot("mask" %in% names(rowData(x)))
  x <- subset(x, !rowData(x)$mask)
  
  message("Checking for NAs (this can take quite a while if HDF5-backed)...") 
  DelayedArray:::set_verbose_block_processing(verbose) # why so slow?
  setAutoBlockSize(1e6) # look at million entries at a time
  t1 <- system.time(naFrac <- rowSums2(is.na(getBeta(x)))/ncol(x))["elapsed"]
  if (verbose) message(sprintf("Computed NA fractions in %.1f seconds", t1))
  rowData(x)$naFrac <- naFrac
  if (all(rowData(x)$naFrac == 0)) {
    if (verbose) message("No unmasked NAs; returning data without imputation.") 
    return(x) 
  }

  # mask excessive NAs 
  maskme <- (rowData(x)$naFrac >= na)
  if (verbose) message("Masking ",sum(maskme)," rows with >= ",na*100,"% NAs.") 
  rowData(x)$mask <- rowData(x)$mask | maskme
  x <- subset(x, !rowData(x)$mask)
  
  # impute any rows where this is practical
  imputable <- rowData(x)$naFrac > 0 & rowData(x)$naFrac < na
  if (verbose) message(sum(imputable), " imputable rows remain.") 
  if (sum(imputable) > 0) {  
    if (sum(imputable) < 10) {
      # grab some nearby probes 
      if (verbose) message("Using additional probes for k-NN imputation.") 
      loci <- subsetByOverlaps(granges(x), resize(granges(x[imputable,]), 5000))
      imputable <- rownames(x) %in% names(loci)
    }
    if (verbose) message("Imputing missing values...") 
    imputeBeta <- function(B) inv.logit(impute.knn(logit(B*sqz))$data)/sqz
    imputed <- imputeBeta(as(getBeta(x[imputable,]), "matrix"))
    assays(x)$Beta[rownames(imputed),] <- imputed
  } else { 
    if (verbose) message("No imputable NAs; returning data without imputation.")
  }

  return(x) 
}
