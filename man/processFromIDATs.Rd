% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/processFromIDATs.R
\name{processFromIDATs}
\alias{processFromIDATs}
\title{read in a directory's worth of IDATs, QC them, optionally sesamize them all}
\usage{
processFromIDATs(frags = 1:3, addgeo = FALSE, justRgSet = FALSE, ...)
}
\arguments{
\item{frags}{which elements of the filenames are relevant? (1:3)}

\item{addgeo}{optional: try to annotate from GEO? (FALSE)}

\item{justRgSet}{optional: dump the rgSet and don't sesamize? (FALSE)}

\item{...}{options to pass to sesame::sesamize}
}
\value{
a GenomicRatioSet (or an rgSet in case of failure)
}
\description{
Note: this function simply reads all IDATs in the current directory,
runs QC, stuffs in metadata, and (if !justRgSet) sesamizes it.
}