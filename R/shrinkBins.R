#' @title Employ an eBayes shrinkage approach for bin-level estimates for A/B inference
#'
#' @description 
#' \code{shrinkBins} returns shrunken bin-level estimates
#'
#' @details This function computes shrunken bin-level estimates using a James-Stein estimator, reformulated as an eBayes procedure
#' 
#' @param x Input SummarizedExperiment object
#' @param original.x Full sample set SummarizedExperiment object
#' @param prior.means The means of the bin-level prior distribution
#' @param chr The chromosome to operate on
#' @param res Resolution to perform the binning
#' @param targets The column/sample/cell names to shrink towards
#' @param assay What assay type this is ("array", "atac", "bisulfite")
#' @param genome What genome are we working with ("hg19", "hg38", "mm9", "mm10")
#'
#' @return A list object to pass to getCorMatrix
#' 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' 
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' #impute to remove NAs
#' imputed.array <- imputeKNN(array.data.chr14, assay = "array")
#' #get the shrunken binned M-values
#' shrunken.bin.array <- shrinkBins(imputed.array, chr = "chr14", assay = "array")
#' 

shrinkBins <- function(x, original.x, prior.means = NULL, chr = NULL,
                       res = 1e6, targets = NULL,
                       assay = c("array", "atac", "bisulfite"),
                       genome = c("hg19", "hg38", "mm9", "mm10")) {
  #match the assay args
  assay <- match.arg(assay)
  
  #match the genome if given
  genome <- match.arg(genome)
  
  #double check the obj class is compatible
  if (!checkAssayType(x)) stop("Input needs to be a SummarizedExperiment")
  
  #get the prior means
  if (is.null(prior.means)) {
    prior.means <- getGlobalMeans(obj=original.x, targets=targets, assay=assay)
  }
  
  #helper function for summary
  atac_fun <- function(x) {
    return(sqrt(mean(x)) * length(x))
  }
  
  #bin the input
  bin.mat <- suppressMessages(switch(assay,
                                     atac = getBinMatrix(x=as.matrix(cbind(assays(original.x)$counts, prior.means)),
                                                         genloc=rowRanges(x), chr=chr, res=res, FUN=atac_fun,
                                                         genome = genome),
                                     array = getBinMatrix(x=as.matrix(cbind(flogit(assays(original.x)$Beta), prior.means)),
                                                         genloc=rowRanges(x), chr=chr, res=res, FUN=median,
                                                         genome = genome),
                                     bisulfite = getBinMatrix(x=as.matrix(cbind(assays(original.x)$counts, prior.means)),
                                                         genloc=rowRanges(x), chr=chr, res=res, FUN=mean,
                                                         genome = genome)))
  
  #shrink the bins using a James-Stein Estimator
  x.shrink <- t(apply(bin.mat$x, 1, function(r) {
    r.samps <- r[!names(r) %in% "globalMean"]
    r.prior.m <- r["globalMean"]
    if (!is.null(targets)) {
      if (length(r.samps[targets]) == 1) {
        stop("Cannot perform targeted bin-level shrinkage with one target sample.")
      }}
    switch(assay,
           atac = .shrinkATAC(x=r.samps, prior=r.prior.m, targets=targets),
           array = .shrinkArrays(x=r.samps, prior=r.prior.m, targets=targets),
           bisulfite = .shrinkBS(x=r.samps, prior=r.prior.m, targets=targets))
    }))
  
  #drop things that are zeroes as global means
  #this can and does crop up in resampling when you have something sparse
  #for instance single-cell data...
  #the correlation will break otherwise
  if (any(bin.mat$x[,"globalMean"] == 0)) {
    bin.mat$gr <- bin.mat$gr[bin.mat$x[,"globalMean"] != 0,]
    x.shrink <- x.shrink[bin.mat$x[,"globalMean"] != 0,]
    bin.mat$x <- bin.mat$x[bin.mat$x[,"globalMean"] != 0,]
  }
  
  return(list(gr=bin.mat$gr, x=x.shrink[,colnames(x)], gmeans=bin.mat$x[,"globalMean"]))
}

#helper functions for computing shrunken means
#shrink bins in (sc)ATAC-seq data
.shrinkATAC <- function(x, prior = NULL, targets = NULL) {
  if (!is.null(targets)) {
    C <- sd(x[targets])
  } else {
    C <- sd(x)
  }
  prior.m <- prior
  #convert back to beta values
  return(prior.m + C*(x - prior.m))
}

#shrink bins in methylation arrays
.shrinkArrays <- function(x, prior = NULL, targets = NULL) {
  if (!is.null(targets)) {
    C <- sd(x[targets])
  } else {
    C <- sd(x)
  }
  prior.m <- prior
  #convert back to beta values
  return(prior.m + C*(x - prior.m))
}

#shrink bisulfite sequencing smoothed M-values
.shrinkBS <- function(x, prior = NULL, targets = NULL) {
  #assumes that M-values exist already
  if (!is.null(targets)) {
    C <- sd(x[targets])
  } else {
    C <- sd(x)
  }
  prior.m <- prior
  return(prior.m + C*(x - prior.m))
}
