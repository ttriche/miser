#' convenience fn
#' 
#' @return  a GRanges
#' 
#' @import Homo.sapiens
#' 
#' @export 
getTxs <- function() { 
  txs <- transcripts(Homo.sapiens, columns=c("SYMBOL","ENSEMBLTRANS"))
  txs <- keepSeqlevels(subset(txs, !duplicated(txs) & !is.na(txs$SYMBOL)),
                       paste0("chr", c(1:22, "X")), pruning="coarse")
  txs$promoter <- resize(flank(txs, 1500), 1700, fix="start")
  return(txs)
}
