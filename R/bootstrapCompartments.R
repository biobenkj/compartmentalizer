#' Non-parametric bootstrapping of compartments and summarization of bootstraps/compute confidence intervals
#'
#' @name bootstrapCompartments
#'
#' @param obj List object of computed compartments for a sample with 'pc' and 'gr' as elements
#' @param original.obj The original, full input SummarizedExperiment of all samples/cells
#' @param bootstrap.samples How many bootstraps to run
#' @param chr Which chromosome to operate on
#' @param assay What sort of assay are we working on
#' @param parallel Whether to run the bootstrapping in parallel
#' @param cores How many cores to use for parallel processing
#' @param targets Targets to shrink towards
#' @param res The compartment resolution
#' @param genome What genome are we working on
#' @param q What sort of confidence intervals are we computing (e.g. 0.95 for 95 percentCI)
#' @param svd The original compartment calls as a GRanges object
#' @param group Whether this is for group-level inference
#' @param bootstrap.means Pre-computed bootstrap means matrix
#'
#' @return Compartment estimates with summarized bootstraps and confidence intervals
#' @import parallel
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 

bootstrapCompartments <- function(obj, original.obj, bootstrap.samples = 1000,
                                  chr = "chr14", assay = c("array", "atac", "bisulfite"),
                                  parallel = TRUE, cores = 2, targets = NULL, res = 1e6,
                                  genome = c("hg19", "hg38", "mm9", "mm10"), q = 0.95,
                                  svd = NULL, group = FALSE, bootstrap.means = NULL) {
  #function for nonparametric bootstrap of compartments and compute 95% CIs
  #check input
  #match the assay args
  assay <- match.arg(assay)

  #double check the obj class is compatible
  if (!checkAssayType(original.obj)) stop("Input needs to be a SummarizedExperiment")
  
  #check the names of the assays
  if (!any(getAssayNames(original.obj) %in% c("Beta", "counts"))) {
    stop("The assay slot should contain either 'Beta' for arrays or 'counts' for atac/bisulfite.")
  }
  
  #if we are using targeted means
  if (!is.null(targets)) original.obj <- original.obj[,targets]
  
  #get the global means we are going to use
  #this could theoretically break if you ask for more bootstraps here than were pre-computed...
  #let's check for one more optimization
  if (bootstrap.samples == ncol(bootstrap.means)) {
    bmeans <- bootstrap.means
  } else {
    bmeans <- sample.int(bootstrap.means, size = bootstrap.samples, replace = FALSE)
    colnames(bmeans) <- rep("globalMean", ncol(bmeans))
    }
  
  #if (ncol(original.obj) < 6) stop("We need more than 5 samples to bootstrap with for the results to be meaningful.")
  if (!parallel) {
    message("Not bootstrapping in parallel will take a long time...")
    #bootstrap and recompute compartments
    resamp.compartments <- lapply(1:ncol(bmeans), function(b) {
      #resample the global means with replacement
      message("Working on bootstrap ", b)
      
      #get the shrunken bins with new global mean
      boot.mean <- as.matrix(bmeans[,b])
      colnames(boot.mean) <- "globalMean"
      s.bins <- shrinkBins(obj, original.obj, prior.means = boot.mean,
                           chr = chr, res = res, assay = assay, genome = genome)
      if (group) cor.bins <- getCorMatrix(s.bins, squeeze = FALSE)
      if (isFALSE(group)) cor.bins <- getCorMatrix(s.bins, squeeze = TRUE)
      #Stupid check for perfect correlation with global mean
      if (any(is.na(cor.bins$binmat.cor))) {
        absig <- matrix(rep(NA, nrow(cor.bins$binmat.cor)))
      } else {
        absig <- getABSignal(cor.bins, assay = assay)
      }
      return(absig)
    })
  } else {
    message("Bootstrapping in parallel with ", cores, " cores.")
    #bootstrap and recompute compartments
    resamp.compartments <- mclapply(1:ncol(bmeans), function(b) {
      #get the shrunken bins with new global mean
      boot.mean <- as.matrix(bmeans[,b])
      colnames(boot.mean) <- "globalMean"
      s.bins <- shrinkBins(obj, original.obj, prior.means = boot.mean,
                           chr = chr, res = res, assay = assay, genome = genome)
      cor.bins <- getCorMatrix(s.bins, squeeze = TRUE)
      #Stupid check for perfect correlation with global mean
      if (any(is.na(cor.bins$binmat.cor))) {
        absig <- matrix(rep(NA, nrow(cor.bins$binmat.cor)))
      } else {
        absig <- getABSignal(cor.bins, assay = assay)
      }
      return(absig)
    }, mc.cores = cores)
  }

  #summarize the bootstraps and compute confidence intervals
  resamp.compartments <- summarizeBootstraps(resamp.compartments, svd,
                                             q = q, assay = assay)
  return(resamp.compartments)
}

#helper function to re-sample
#this was inspired by https://github.com/sgibb/bootstrap/blob/master/R/helper-functions.R
.resampleMatrix <- function(x, size=ncol(x)) {
  samp.to.select <- sample.int(ncol(x), size=size, replace=TRUE)
  return(x[, samp.to.select])
}
