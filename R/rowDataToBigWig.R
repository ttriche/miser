#' convenience function to dump (e.g.) W matrix values to bigWig files
#'
#' @param x    a GRanges, a SummarizedExperiment-like object, or its rowRanges
#' @param mcol  the name of the mcols column in x or rowData(x) to bigWigify
#' @param fname the name of the bigWig file (else default to mcol.genome.bw) 
#' @param seqs  optional seqinfo for any `x` without seqlengths
#' @param nmf   NMF object with rownames identical to x, for its W matrix
#' 
#' @import rtracklayer
#' @import RcppML
#'
#' @export
rowDataToBigWig <- function(x, mcol, fname=NULL, seqs=NULL, nmf=NULL) { 

  if (is.null(fname)) fname <- paste(mcol, unique(genome(x)), "bw", sep=".")
  if (!is.null(seqs)) seqlengths(x) <- seqlengths(seqs)[seqlevels(x)]
  rtracklayer::export(.prepGRanges(x, mcol=mcol, nmf=nmf), fname)

}


# helper fn
.prepGRanges <- function(gr, mcol, nmf=NULL) {

  gr <- .nmfAsGRanges(granges(gr), mcol=mcol, nmf=nmf)
  stopifnot(mcol %in% names(mcols(gr)))
  score(gr) <- mcols(gr)[, mcol]
  return(gr)

}


# helper fn
.nmfAsGRanges <- function(gr, mcol=NULL, nmf=NULL) { 

  if (!is.null(nmf)) {
    stopifnot(all(rownames(nmf@w) %in% names(gr))) 
    rows <- intersect(rownames(nmf@w), names(gr))
    gr <- gr[rows]
    if (!is.null(mcol)) {
      mcols(gr) <- nmf@w[rows, mcol]
    } else { 
      mcols(gr) <- nmf@w[rows, ]
    }
  }
  return(gr)

}
