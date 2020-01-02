#' Dirichlet smoothing of BS-seq counts prior to compartment inference
#'
#' @name smoothBSCounts
#'
#' @description helper function for compartment inference (shrink-by-smoothing logit frac5mC)
#' 
#' We want something with nominally Gaussian error for compartment inference, so
#' this function grabs suitable (default >= 3 reads in >=2 sample) measurements
#' and turns them into lightly moderated, logit-transformed methylated-fraction
#' estimates (also known, unfortunately, as M-values) for compartment calling,
#' by performing Dirichlet smoothing (adding `k` reads to M and U support).
#' 
#' @param obj         a BSseq object with methylated and total reads 
#' @param minCov      minimum read coverage for landmarking samples (3)
#' @param minSamp     minimum landmark samples with >= minCov (2)
#' @param k           pseudoreads for smoothing (0.5)
#' @param r           regions to collapse over (default is NULL, do it by CpG)
#' 
#' @return            smoothed logit(M/Cov) matrix with coordinates as row names
#' 
#' @import            gtools
#' @import            bsseq
#'
#' @export
#' 
#' @examples 
#' 

smoothBSCounts <- function(obj, minCov=3, minSamp=2, k=0.5, r=NULL) {
  #check if the input is sane
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")
  
  #make sure M and Cov exist
  if (!any(getAssayNames(obj) %in% c("M", "Cov"))) {
    stop("The assays slot of the bsseq object needs to have 'M' and 'Cov'.")
  }
  return(.getLogitFracMeth(obj, minCov = minCov, minSamp = minSamp,
                           k = k, r = r))
}

#helper function to calculate smoothed counts
.getLogitFracMeth <- function(x, minCov=3, minSamp=2, k=0.1, r=NULL) {
  
  # do any loci/regions have enough read coverage in enough samples? 
  if (!is.null(r) && is(r, "GenomicRanges")) {
    covgs <- getCoverage(x, sort(r), type="Cov", what="perRegionTotal")
  } else { 
    covgs <- getCoverage(x, type="Cov", what="perBase")
  } 
  
  usable <- DelayedMatrixStats::rowSums2(covgs >= minCov) >= minSamp
  if (!any(usable)) stop("No usable loci/regions ( >= minCov in >= minSamp )!")
  
  # construct a subset of the overall BSseq object with smoothed mvalues 
  if (!is.null(r) && is(r, "GenomicRanges")) {
    smoothed.counts <- .getSmoothedLogitFrac(x, k=k, minCov=minCov, r=subset(sort(r), usable))
    return(assays(x)$counts <- smoothed.counts)
  } else { 
    smoothed.counts <- .getSmoothedLogitFrac(subset(x, usable), k=k, minCov=minCov)
    #convert to GRanges
    smooth.ranges <- granges(as.matrix(smoothed.counts))
    #make into a SummarizedExperiment object
    smooth.se <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(smoothed.counts)),
                                      rowRanges=granges(smooth.ranges),
                                      colData=colData(x))
    return(smooth.se)
  } 
  
}

# helper fn
.getSmoothedLogitFrac <- function(x, k=0.1, minCov=3, maxFrac=0.5, r=NULL) {
  
  if (!is.null(r) && is(r, "GenomicRanges")) {
    M <- getCoverage(x, sort(r), type="M", what="perRegionTotal")
    U <- getCoverage(x, sort(r), type="Cov", what="perRegionTotal") - M 
    rnames <- as.character(sort(r))
  } else { 
    M <- getCoverage(x, type="M", what="perBase")
    U <- getCoverage(x, type="Cov", what="perBase") - M 
    rnames <- as.character(granges(x))
  } 
  
  res <- logit((M + k) / ((M + k) + (U + k))) 
  rownames(res) <- rnames 
  
  makeNA <- ((M + U) < minCov)
  maxPct <- paste0(100 * maxFrac, "%")
  tooManyNAs <- (DelayedMatrixStats::colSums2(makeNA)/nrow(x)) > maxFrac
  if (any(tooManyNAs)) {
    message(paste(colnames(x)[tooManyNAs],collapse=", ")," are >",maxPct," NA!")
  }
  res[ makeNA ] <- NA
  return(res)
  
}
