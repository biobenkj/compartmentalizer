% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shrinkBins.R
\name{shrinkBins}
\alias{shrinkBins}
\title{Employ an eBayes shrinkage approach for bin-level estimates for A/B inference}
\usage{
shrinkBins(
  x,
  original.x,
  prior.means = NULL,
  chr = NULL,
  res = 1000000,
  targets = NULL,
  jse = FALSE,
  assay = c("array", "atac", "bisulfite", "rna"),
  genome = c("hg19", "hg38", "mm9", "mm10")
)
}
\arguments{
\item{x}{Input SummarizedExperiment object}

\item{original.x}{Full sample set SummarizedExperiment object}

\item{prior.means}{The means of the bin-level prior distribution}

\item{chr}{The chromosome to operate on}

\item{res}{Resolution to perform the binning}

\item{targets}{The column/sample/cell names to shrink towards}

\item{jse}{Whether to use a James-Stein estimator instead}

\item{assay}{What assay type this is ("array", "atac", "bisulfite", "rna")}

\item{genome}{What genome are we working with ("hg19", "hg38", "mm9", "mm10")}
}
\value{
A list object to pass to getCorMatrix
}
\description{
\code{shrinkBins} returns shrunken bin-level estimates
}
\details{
This function computes shrunken bin-level estimates using a James-Stein estimator, reformulated as an eBayes procedure
}
\examples{
data("meth_array_450k_chr14", package = "compartmap")
#impute to remove NAs
imputed.array <- imputeKNN(array.data.chr14, assay = "array")
#get the shrunken binned M-values
shrunken.bin.array <- shrinkBins(imputed.array, chr = "chr14", assay = "array")

}
