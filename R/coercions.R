
# need to work out namespace/sealing issues here
if (FALSE) { 
  
  
setAs("RGChannelSet", "SummarizedExperiment", 
      function(from) {
        metadata(from)$annotation <- annotation(from)
        SummarizedExperiment(assays = assays(from),
                             colData = colData(from),
                             rowData = rowData(from), 
                             metadata = metadata(from))
      })


setAs("SummarizedExperiment", "RGChannelSet", 
      function(from) {
        stopifnot(all(c("Green", "Red") %in% assayNames(from)))
        stopifnot("annotation" %in% names(metadata(from)))
        RGChannelSet(Green = getGreen(from),
                     Red = getRed(from),
                     annotation = metadata(from)$annotation,
                     colData = colData(from),
                     rowData = rowData(from), 
                     metadata = metadata(from))
      })


setAs("GenomicRatioSet", "SummarizedExperiment", 
      function(from) {
        metadata(from)$annotation <- annotation(from)
        metadata(from)$preprocessMethod <- preprocessMethod(from)
        SummarizedExperiment(assays = assays(from),
                             colData = colData(from),
                             rowData = rowData(from), 
                             metadata = metadata(from))
      })


setAs("SummarizedExperiment", "GenomicRatioSet", 
      function(from) {
        stopifnot(any(c("Beta", "M") %in% assayNames(from)))
        stopifnot("preprocessMethod" %in% names(metadata(from)))
        stopifnot("annotation" %in% names(metadata(from)))
        GenomicRatioSet(Beta = assays(from, "Beta"), 
                        CN = assays(from, "CN"), 
                        M = assays(from, "M"),
                        colData = colData(from),
                        rowData = rowData(from), 
                        metadata = metadata(from),
                        annotation = metadata(from)$annotation,
                        preprocessMethod = metadata(from)$preprocessMethod)
      }) 


} 
