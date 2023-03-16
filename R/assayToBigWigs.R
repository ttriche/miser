#' Dump data as bigWig tracks for (e.g.) igv.js (which can't handle multiGCTs)
#'
#' @param x         a RangedSummarizedExperiment-like object of some sort
#' @param assay     name or index of an assay to dump (default: assays(x)[[1]])
#' @param genome    (optional) genome to use (default: look in genome(x))
#' @param onlycg    (optional) only dump CpG probes (for methylation)? (TRUE)
#'
#' @return list     a list of filenames dumped
#'
#' @import rtracklayer 
#' 
#' @export
assayToBigWigs <- function(x, assay=1, genome=NULL, onlycg=TRUE) { 

  if (is.null(genome)) genome <- unique(genome(x))
  if (length(genome) == 0) stop("No genome specified or provided")
  if (length(genome) > 1) stop("Multiple genomes detected or provided")
  if (is.null(assay)) assay <- 1 
  if (onlycg & length(grep("^cg", rownames(x))) > 1) {
    x <- x[grep("^cg", rownames(x)), ] 
  }

  gr <- granges(x)
  for (i in colnames(x)) { 
    gr$score <- getBeta(x)[, i]
    filename <- paste(colnames(x)[i], g, "bw", sep=".")
    export(gr, filename)
    filenames <- append(filenames, filename)
  }

  invisible(filenames)

}
