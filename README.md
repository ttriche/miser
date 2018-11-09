# miser

[![Build Status](https://travis-ci.org/ttriche/miser.png?branch=master)](https://travis-ci.org/ttriche/miser)  [![codecov](https://codecov.io/gh/ttriche/miser/branch/master/graph/badge.svg)](https://codecov.io/gh/ttriche/miser)

## Batch-preprocess data from GEO and other sources

```r

library(miser)

# in a clean new directory...
GSMs <- paste0("GSM230915", 4:7)
IDATs <- getGSMs(GSMs, cachePath=".") # downloads metadata here 
grSet <- sesamizeGEO(IDATs, annot=TRUE, HDF5=TRUE, cachePath=".")
colData(grSet) # autopopulated from GEOmetadb when annot=TRUE
```

## Save the (joint or individual) data as HDF5 files (out-of-core computation)

```r

```

## Call differentially methylated regions out-of-core using HDF5 backing

```r

```

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
