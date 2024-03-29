% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dumpQCfiles.R
\name{dumpQCfiles}
\alias{dumpQCfiles}
\title{dump QC (NAs, flags, XY stats) and (optionally) SNP calls and beta values}
\usage{
dumpQCfiles(grSet, stub = NULL, path = ".", snps = TRUE, betas = FALSE)
}
\arguments{
\item{grSet}{GenomicRatioSet with betas and QC metrics to dump out}

\item{stub}{mandatory prefix for derived files (e.g. stub_QC.csv.gz)}

\item{path}{where to save the dumped CSV file\link{s} (".")}

\item{snps}{dump SNPs and SNP calls? (TRUE)}

\item{betas}{dump beta values too? (FALSE)}
}
\value{
\if{html}{\out{<div class="sourceCode">}}\preformatted{       the status of the final operation
}\if{html}{\out{</div>}}
}
\description{
NA fraction, median CN, control flags, XY stats are saved in stub_QC.csv.gz
SNP calls are saved in stub_SNPcalls.csv.gz if snps == TRUE
SNP intensities are saved in stub_SNPs.csv.gz if snps == TRUE
Sesamized betas are saved in stub_betas.csv.gz if betas == TRUE
}
\seealso{
\if{html}{\out{<div class="sourceCode">}}\preformatted{      SNPcalls
}\if{html}{\out{</div>}}
}
