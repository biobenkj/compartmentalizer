#' @title Filter and impute missing values prior to A/B compartment inference with k-nearest neighbor
#' 
#' @description 
#' \code{filterAndImputeWGBS} returns filtered and smoothed signal from whole genome bisulfite sequencing data
#'
#' @details 
#' This function filters and performs Dirichlet smoothing of WGBS signal prior to A/B compartment inference
#'
#' @param obj Input BSseq object
#' @param k The pseudocount to be added for smoothing - default uses Jeffery's prior
#' @param minCov Minimum coverage of a locus to be smoothed
#' @param minSamp Minimum samples that have sufficient coverage to be smoothed
#' @param rowmax The maximum percent missing data allowed per loci - default = 0.5
#' @param colmax The maximum percent missing data allowed per sample - default = 0.8
#' @param openseafilt Whether to filter loci to just open sea CpG loci
#' @param genome The genome to analyze - currently only supports hg19
#'
#' @return Filtered and imputed matrix 
#' @import biscuiteer
#' @import impute
#' @import Homo.sapiens
#'
#' @export
#' 
#' @examples 
#' 
#' library(biscuiteer)
#' library(Homo.sapiens)
#' library(impute)
#' 
#' #Load in some example data
#' data(cell_cycle_hansen_chr22, package = "compartmap")
#' 
#' #Filter and impute missing data
#' filt.imputed.chr22 <- filterAndImputeWGBS(data.chr22)

#Filter and impute missing values prior to binning
filterAndImputeWGBS <- function(obj, k = 0.5, minCov = 3, minSamp = 2, rowmax = 0.5, colmax = 0.8, openseafilt = FALSE, genome = "hg19") {
  #Check whether we should filter to opensea CpG loci
  #Currently only supports hg19
  if (openseafilt) obj <- .filtOpenSeaLoci(obj, genome)
  #Coerce BSseq object to matrix
  if (is(obj, "BSseq")) obj <- .coerceBSseq(obj, k = k, minCov = minCov, minSamp = minSamp)
  message(paste0("Removing loci with greater than ", rowmax*100, "% NA values prior to knn imputation..."))
  obj.clean <- obj[rowMeans(is.na(obj)) < rowmax,]
  message("Starting knn imputation...")
  #This will spew all kinds of output but some users like seeing things are working I suppose...
  obj.clean.imputed <- impute.knn(obj.clean, rowmax = rowmax, colmax = colmax)$data
  return(obj.clean.imputed)
}

#Helper function to subset the input BSseq object to open sea CpG probes
#Will eventually support hg38 and other genomes 
.filtOpenSeaLoci <- function(obj, genome = "hg19") {
  #Load the openSeas.hg19 data
  data("openSeas.hg19", package = "biscuiteer")
  #Subset by overlaps
  message("Filtering to open sea CpG loci...")
  obj.openseas <- subsetByOverlaps(obj, openSeas.hg19)
  return(obj.openseas)
}

#Helper function to smooth WGBS data prior to imputation
#Realize and coerce to a matrix from bsseq object
.coerceBSseq <- function(bsseq, k = 0.5, minCov = 3, minSamp = 2) {
  if (is(bsseq, "BSseq")) {
    message("Assuming this hasn't been smoothed yet since it is still a BSseq object...")
    message("Smoothing...")
    bsseq.smooth <- getLogitFracMeth(bsseq, k = k, minCov = minCov, minSamp = minSamp)
    message("Realizing DelayedArray object...")
    message("This could take a little while...")
    bsseq.smooth <- realize(bsseq.smooth)
    message("Coercing to a matrix for compartment inference...")
    bsseq.obj <- as.matrix(bsseq.smooth)
  }
  else if (is(bsseq, "matrix")) {
    message("Realizing DelayedArray object...")
    message("This could take a little while...")
    #Add a try statement here in case it's just a standard matrix at this point
    bsseq <- try(realize(bsseq))
    message("Trying to coerce to a matrix for compartment inference...")
    bsseq.obj <- as.matrix(bsseq)
  }
  else (stop(paste0("Not sure what do with ", class(bsseq))))
  return(bsseq.obj)
}
