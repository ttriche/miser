#' Dump DNA methylation values for quick and dirty IGV viewing
#'
#' @param x         a SummarizedExperiment-like object of some sort
#' @param assay     the name of the assay to dump (default: assays(x)[[1]])
#' @param group     a grouping factor to dump by (default: dump everyone) 
#' @param stub      the filename stub to dump (default: name of x, or x_group)
#' @param genome    the genome to use (default is hg19, due to annotations)
#' 
#' @return list     a list of filenames dumped
#'
#' @export
asGCT <- function(x, assay=1, group=NULL, stub=NULL, genome="hg19") { 

  x <- sort(x) ## must be sorted by position for GCT and toTDF 
  if (is.null(stub)) stub <- as.character(match.call()["x"])
  if (is.null(assay)) assay <- 1 

  gr <- granges(x)
  if (is.null(group)) { 
    ## dump everyone at the same time 
    filename <- .dumpGCT(assays(x)[[assay]], gr, stub, genome)
    filenames <- list(stub=filename)
  } else { 
    ## dump by group
    filenames <- list()
    for (g in levels(as.factor(group))) {
      grp <- which(group == g)
      grpstub <- paste(stub, g, sep=".")
      filename <- .dumpGCT(assays(x[, grp])[[assay]], gr, grpstub, genome)
      filenames <- append(filenames, filename)
    }
  }

  invisible(filenames)

}


# add coordinates for features 
.dumpGCT <- function(mat, gr, stub, genome="hg19") {
  fname <- paste(stub, genome, "gct", sep=".")
  description <- paste0("na |@", as.character(gr), "|")
  contents <- data.frame(NAME=rownames(mat), Description=description, mat)
  message("Dumping values to ", fname, "...")
  ## formatting for GCT 
  cat("#1.2", "\n", file=fname)
  cat(nrow(mat), "\t", ncol(mat), "\n", file=fname, append=TRUE) 
  cat(colnames(mat), "\n", file=fname, sep="\t", append=TRUE) 
  write.table(contents,
              file=fname, append=T, row.names=F, col.names=F, quote=F, sep="\t")
  return(fname)
}


# as above
.toBED4 <- function(gr, filename) {
  coord <- as.data.frame(gr)[,1:3]
  coord[,2] <- coord[,2] - 1 ## derrrrp
  coord <- cbind(coord, rownames(coord))
  write.table(coord, row.names=F, col.names=F, quote=F, file=filename, sep="\t")
}
