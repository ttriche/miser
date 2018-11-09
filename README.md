# miser

[![Build Status](https://travis-ci.org/ttriche/miser.png?branch=master)](https://travis-ci.org/ttriche/miser)  [![codecov](https://codecov.io/gh/ttriche/miser/branch/master/graph/badge.svg)](https://codecov.io/gh/ttriche/miser)

## Batch-preprocess data from GEO and other sources

```r
library(miser)

# in a clean new directory...
GSMs <- paste0("GSM230915", 4:7)
IDATs <- getGSMs(GSMs, cachePath=".") # downloads metadata here 
grSet <- sesamizeGEO(IDATs, annot=TRUE, HDF5=TRUE, cachePath=".")

# Since annot=TRUE, the column data for grSet is automatically filled in.
# It can take a while to download the GEO metadatabase required for this.
colData(grSet)
```

## Save the result as HDF5 files (for out-of-core computation)

This GenomicRatioSet is tiny, but once you start manipulating thousands of 
samples at a time (e.g. the NOPHO ALL data or DKFZ nervous system tumors), 
it comes in handy to only load data into RAM as necessary. Using HDF5 files
to hold the data on disk allows R not to be such a pig about memory, in 
exchange for somewhat slower computation (since disk or cache needs access).  

```r
grSet <- saveAsHDF5(grSet, "GSM230915") # now HDF5-backed from GSM230915 folder
```

By assigning the result of saveAsHDF5 to grSet, we are switching to disk-backed
out-of-core storage. This is usually a good thing, but can be slow if package 
authors aren't aware of how access works. We rewrote a simple DMRcate wrapper
as an example of how a common task can be adapted to out-of-core computation.


## Call differentially methylated regions out-of-core using HDF5 backing

In this example we have primary prostate tissue and a prostate cancer cell line,
so we evaluate differentially methylated regions between replicates of each. 
Missing or masked probes are imputed using k-nearest-neighbors (kNN) by default,
and the X and Y chromosome would usually be dropped (in this case, since both 
the cell line and the control come from a male, it might be better to keep XY). 

```r
design <- with(colData(grSet), model.matrix( ~ tissue))
results <- getDMRs(grSet, design) 
results 
results$DMRs 
bySign <- split(results$DMRs, sign(results$DMRs$meanbetafc))
sapply(bySign, length)
sapply(sapply(bySign, width), summary)
```

It is somewhat unusual to see shorter regions hypomethylated than 
hypermethylated, when comparing a cancer cell line to primary tissue,
but that is what the summaries by sign indicate. It would be worth 
annotating these to see where they fall in terms of chromatin states,
as with `erma` (example to be added). 

## A/B compartments

The `compartmap` package (part of Bioconductor) can be used to infer chromatin
compartments (open/closed regions) from DNA methylation microarrays. As far as 
we know, it works pretty much the same with HDF5 or in-core storage. 

## Other stuff 

Gene set analysis, enrichment tests, and so forth can be done in `missMethyl`.

## Installing

```r
install.packages("BiocManager") 
BiocManager::install("ttriche/miser") 
```

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
