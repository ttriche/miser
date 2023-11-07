#' "fix" a SummarizedExperiment (for which IDATs may be unavailable) with Sesame
#' 
#' @param x         SummarizedExperiment-derived object with assays Red & Green
#' @param genome    Which genome should the probes be mapped to? (hg19)
#' @param BPPARAM   a BiocParallelParam object (SerialParam())
#' @param mft       a data.frame with columns Probe_ID, M, U, and col (guess)
#' @param sce       return SingleCellExperiment rather than a RatioSet? (TRUE)  
#' @param verbose   be verbose? (FALSE) 
#' @param ...       additional arguments passed to openSesame
#'
#' @return a sesamized SingleCellExperiment or GenomicRatioSet from `x`
#'
#' @import BiocParallel
#' @import SummarizedExperiment 
#' @import SingleCellExperiment 
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @importFrom BiocManager available 
#' 
#'
#' @examples
#' \donttest{
#'
#'   # reprocess some 450k data from an existing RGChannelSet 
#'   if (require(FlowSorted.CordBloodNorway.450k) && 
#'     require(IlluminaHumanMethylation450kmanifest) && 
#'     require(IlluminaHumanMethylation450kanno.ilmn12.hg19)) {
#'     sesamized <- sesamize(FlowSorted.CordBloodNorway.450k[,1:2])
#'   } 
#'
#'   # reprocess some TCGA 450k data from IDATs 
#'   if (require(minfi) && 
#'     require(TCGAMethylation450k) && 
#'     require(IlluminaHumanMethylation450kmanifest) && 
#'     require(IlluminaHumanMethylation450kanno.ilmn12.hg19)) {
#'     TCGA_path <- system.file("extdata","idat", package="TCGAMethylation450k")
#'     TCGA_IDATs <- list.files(TCGA_path)
#'     TCGA_stubs <- unique(sub("_(Red|Grn).idat", "", TCGA_IDATs))
#'     TCGA_targets <- data.frame(Basename=TCGA_stubs)
#'     TCGA_rgSet <- read.metharray.exp(base=TCGA_path, targets=TCGA_targets)
#'     TCGA_grSet <- sesamize(TCGA_rgSet)
#'   }
#'
#' }
#'
#' @export 
#'
sesamize <- function(x, genome="hg19", BPPARAM=SerialParam(), mft=NULL, sce=TRUE, verbose=FALSE, ...) { 

  stopifnot(is(x, "SummarizedExperiment"))
  stopifnot(all(c("Green", "Red") %in% assayNames(x)))

  if (ncol(x) > 1) {

    # {{{ recursively sesamize each sample
    nms <- colnames(x)
    names(nms) <- nms
    res <- do.call(cbind, 
                   bplapply(nms, 
                            function(y) sesamize(x[,y], mft=mft), 
                            BPPARAM=BPPARAM))
    ppm <- c(rg.norm="SeSAMe (type I)",
             p.value="SeSAMe (pOOBAH)",
             manifest=annotation(x)["array"])

    if (!sce) {
      # {{{ deprecated 
      res <- minfi::mapToGenome(res)
      mfst <- packageVersion(paste(minfi::annotation(x), collapse="anno."))
      res@preprocessMethod <- ppm
      # }}} 
    }

    metadata(res) <- metadata(x)
    metadata(res)$annotation <- annotation(x)
    metadata(res)$preprocessMethod <- ppm
    colData(res) <- colData(x)
    res <- res[which(rowSums(is.na(getBeta(res))) < ncol(res)), ]
    res <- res[, which(colSums(is.na(getBeta(res))) < nrow(res))]
    if (sce) {
      res <- as(res, "SingleCellExperiment")
      rowRanges(res) <- .getMappings(x)[rownames(res)]
    }
    return(res)

  } else {
   
    # sesamize one sample. this is the crux  
    if (verbose) message("Sesamizing ", colnames(x), "...")
    rg <- cbind(G=getGreen(x), R=getRed(x)) # newly generic 
    colnames(rg) <- c("G", "R")
    if (is.null(mft)) { # kludge; ask Wanding how to do better
      data("mfts", package="miser")
      if (nrow(rg) > 1000000) {
        mft <- mfts$EPIC$ordering
      } else if (nrow(rg) > 600000) {
        mft <- mfts$HM450$ordering
      } else if (nrow(rg) < 600000) {
        mft <- mfts$MM285$ordering
      }
    }
    attr(rg, "platform") <- sub("HMEPIC", "EPIC", 
                                sub("IlluminaHumanMethylation", "HM", 
                                    sub("k$", "", annotation(x)["array"])))
    sdf <- sesame:::chipAddressToSignal(rg, mft=mft) # see above for kludge
    if (verbose) message(nrow(sdf), " probe addresses for ", colnames(x))
    Beta <- matrix(openSesame(sdf), ncol=1, 
                   dimnames=list(sdf$Probe_ID, colnames(x)))
    CN <- matrix(log2(totalIntensities(sdf))[rownames(Beta)], ncol=1, 
                 dimnames=list(sdf$Probe_ID, colnames(x)))
    res <- SummarizedExperiment(assays=list(Beta=Beta, CN=CN)) 
    if (verbose) message(nrow(res), " probes processed for ", colnames(x))
    preprocessMethod(res) <- preprocessMethod(x)
    annotation(res) <- annotation(x)
    return(res)

  }
}


# helper fn
.getAry <- function(ary) {
  if (is(ary, "SummarizedExperiment")) ary <- annotation(ary)["array"]
  sub("450", "HM450", sub("IlluminaHumanMethylation", "", sub("k$", "", ary)))
}


# helper fn
.getPossible <- function(ary) {
  sesameDataList(filter=paste(.getAry(ary), "probeInfo", sep="."))$Title
}


# helper fn
.getMappings <- function(ary, genome="hg19") { 
  possible <- .getPossible(ary)
  if (length(possible) > 1) {
    message("Multiple possible mappings:")
    for (i in possible) message(i)
    stop("Please choose one manually.")
  } else { 
    sdg <- sesameDataGet(possible)
    mapping <- paste("mapped", "probes", genome, sep=".")
    mapped <- sdg[[mapping]]
    genome(mapped) <- genome
    
    # completely ignoring the existence of MM285...
    masking <- ifelse(genome == "hg38", "mask", paste("mask", genome, sep="."))
    masked <- sdg[[masking]]
    mcols(mapped)[["mask"]] <- names(mapped) %in% masked
  }
  return(mapped)
}
