#' run (e.g.) eBayes(lmFit(assay(x), design)) on columns of (e.g.) assay(x)
#' 
#' if `y` is specified, use that instead of assay(x) but constrain it the same
#' if `purify` is specified, purify assay(x) before proceeding with the fit(s)
#'
#' @param   x         a SummarizedExperiment-like object
#' @param   columns   colData(x) column names to regress as predictors
#' @param   y         optional matrix w/ ncols(x) columns (default is assay(x))
#' @param   how       "lmFit" (default) or "concordance" (survival::concordance)
#' @param   check     ensure that the selected columns of y are non-NA? (TRUE) 
#' @param   purify    purify `x`? (default is FALSE, but this can be useful) 
#' @param   purity    if `purify` == TRUE, a required vector of purities for x
#' @param   betas     if `purify` == TRUE, are the factors beta values? (FALSE)
#' @param   ...       additional arguments to lmFit, if how == "fit"
#'
#' @return            a fitted model (use topTable on it) or concordance object
#'
#' @seealso           purify
#' @seealso           limma::lmFit
#' @seealso           limma::topTable
#' @seealso           survival::concordance
#'
#' @import            survival
#' @import            limma
#'
#' @export
regressFactors <- function(x, columns, y=NULL, how=c("lmFit","concordance"), check=TRUE, purify=FALSE, purity=NULL, betas=FALSE, ...) {

  stopifnot(all(columns %in% names(colData(x))))
  if (!is.null(y)) stopifnot(ncol(y) == ncol(x))
  if (is.null(y)) y <- assay(x)
  if (purify) y <- purify(y, purity=purity, betas=betas)

  preds <- colData(x)[, columns]
  keep <- which(rowSums(is.na(preds)) == 0)
  message("Kept ", length(keep), " out of ", ncol(x), " samples.")
  stopifnot(length(keep) >= 3)
  preds <- preds[keep, ]
  
  y <- y[, keep] 
  if (sum(is.na(y)) > 0) stop("NA values present in y; cannot fit.")
  
  how <- match.arg(how)
  if (how == "concordance") stop("`concordance` is not yet supported.") 

  design <- with(preds, model.matrix(~ .))
  message("coefs (for topTable(fit)): ")
  for (i in seq_along(colnames(design))) message(i, ": ", colnames(design)[i])
  eBayes(lmFit(y, design))

}
