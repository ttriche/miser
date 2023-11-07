# Make generic SummarizedExperiment-derived objects behave as if minfi objects.

if (!isGeneric("getGreen")) setGeneric("getGreen")
setMethod("getGreen", "SummarizedExperiment", 
          function(object) assay(object, "Green"))

if (!isGeneric("getRed")) setGeneric("getRed")
setMethod("getRed", "SummarizedExperiment", 
          function(object) assay(object, "Red"))

if (!isGeneric("annotation")) setGeneric("annotation")
setMethod("annotation", "SummarizedExperiment", 
          function(object, ...) metadata(object)$annotation)

if (!isGeneric("annotation<-")) setGeneric("annotation<-")
setReplaceMethod("annotation", 
                 c(object="SummarizedExperiment", value="ANY"),
                 function(object, ..., value) {
                   metadata(object)$annotation <- value
                   return(object)
                 })

if (!isGeneric("preprocessMethod")) setGeneric("preprocessMethod")
setMethod("preprocessMethod", "SummarizedExperiment", 
          function(object) metadata(object)$preprocessMethod)

if (!isGeneric("preprocessMethod<-")) {
  setGeneric("preprocessMethod<-",
             function(object, ..., value) standardGeneric("preprocessMethod<-"))
}

setReplaceMethod("preprocessMethod", 
                 c(object="SummarizedExperiment", value="ANY"),
                 function(object, ..., value) {
                   metadata(object)$preprocessMethod <- value
                   return(object)
                 })

# helper fn
.ppm <- function(object, ..., value) {
  object@preprocessMethod <- value
  return(object)
}
setReplaceMethod("preprocessMethod", c(object="MethylSet", value="ANY"), .ppm)
setReplaceMethod("preprocessMethod", c(object="RatioSet", value="ANY"), .ppm)

# helper fn
.M2B <- function(x) (2 ** x) / (1 + (2 ** x))

# helper fn
.B2M <- function(p) log2(p / (1 - p))

# using the above 
if (!isGeneric("getBeta")) setGeneric("getBeta")
setMethod("getBeta", "SummarizedExperiment", 
          function(object, ...) {
            if ("Beta" %in% assayNames(object)) assay(object, "Beta")
            else if ("M" %in% assayNames(object)) .M2B(assay(object, "M"))
            else stop("No Beta or M assay found")
          })

# also using the above
if (!isGeneric("getM")) setGeneric("getM")
setMethod("getM", "SummarizedExperiment", 
          function(object, ...) {
            if ("M" %in% assayNames(object)) assay(object, "M")
            else if ("Beta" %in% assayNames(object)) .B2M(assay(object, "Beta"))
            else stop("No Beta or M assay found")
          })

if (!isGeneric("getCN")) setGeneric("getCN")
setMethod("getCN", "SummarizedExperiment", 
          function(object, ...) assay(object, "CN"))

if (!isGeneric("getMeth")) setGeneric("getMeth")
setMethod("getMeth", "SummarizedExperiment", 
          function(object) getBeta(object) * (2 ** getCN(object)))

if (!isGeneric("getUnmeth")) setGeneric("getUnmeth")
setMethod("getUnmeth", "SummarizedExperiment", 
          function(object) (1 - getBeta(object)) * (2 ** getCN(object)))
