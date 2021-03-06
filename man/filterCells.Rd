% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterCells.R
\name{filterCells}
\alias{filterCells}
\title{Filter cells using either k-means or Dirichlet process means clustering of sparsity metrics}
\usage{
filterCells(
  sparsity.mat,
  rse.obj,
  cluster.method = c("kmeans", "dpmeans"),
  clusters = NULL,
  tol = 0.1,
  plot.data = FALSE,
  invert = FALSE
)
}
\arguments{
\item{sparsity.mat}{A matrix of summarized sparsity measures}

\item{rse.obj}{The unfiltered RangedSummarizedExperiment object}

\item{cluster.method}{Clustering method to use (default: kmeans)}

\item{clusters}{How many clusters to generate; if NULL it will autopick the cluster number (default: NULL)}

\item{tol}{The tolerance or minimum difference in fraction of between cluster sum of squares over total for k-means auto-picking cluster number (default: 0.1)}

\item{plot.data}{Whether to plot the data}

\item{invert}{Invert which cluster is used to filter}
}
\value{
A filtered RangedSummarizedExperiment object and/or plot of the filtered data
}
\description{
Filter cells using either k-means or Dirichlet process means clustering of sparsity metrics
}
\examples{

# FIXME: add example

}
