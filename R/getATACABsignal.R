#' @title Estimate A/B compartments from ATAC-seq data
#'
#' @description 
#' \code{getATACABsignal} returns estimated A/B compartments from ATAC-seq data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param parallel Whether to run samples in parallel
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param targets Samples/cells to shrink towards
#' @param cores How many cores to use when running samples in parallel
#' @param bootstrap Whether we should perform bootstrapping of inferred compartments
#' @param num.bootstraps How many bootstraps to run
#' @param genome What genome to work on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#' @param group Whether to treat this as a group set of samples
#' @param boot.parallel Whether to run the bootstrapping in parallel
#' @param boot.cores How many cores to use for the bootstrapping
#'
#' @return A RaggedExperiment of inferred compartments
#' @import SummarizedExperiment
#' @import parallel
#' @import RaggedExperiment
#' @import Homo.sapiens
#' @import Mus.musculus
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @import BSgenome.Mmusculus.UCSC.mm9
#' @export
#' @examples
#' data("groupATAC_raw_filtered_chr14", package = "compartmap")
#' atac_compartments <- getATACABsignal(filtered.data.chr14, parallel=F, chr="chr14", bootstrap=F, genome="hg19")

getATACABsignal <- function(obj, res = 1e6, parallel = TRUE, chr = NULL,
                             targets = NULL, cores = 2,
                             bootstrap = TRUE, num.bootstraps = 1000,
                             genome = c("hg19", "hg38", "mm9", "mm10"),
                             other = NULL, group = FALSE,
                             boot.parallel = TRUE, boot.cores = 2) {
  
  #gather the chromosomes we are working on
  if (is.null(chr)) {
    message("Assuming we want to process all chromosomes.")
    #get what chromosomes we want
    chr <- getChrs(obj)
  }
  
  #get the column names
  if (is.null(colnames(obj))) stop("colnames needs to be sample names.")
  columns <- colnames(obj)
  names(columns) <- columns
  
  #precompute global means
  prior.means <- getGlobalMeans(obj = obj, targets = targets, assay = "atac")
  
  if (bootstrap) {
    message("Pre-computing the bootstrap global means.")
    bmeans <- precomputeBootstrapMeans(obj = obj, targets = targets, num.bootstraps = num.bootstraps,
                                       assay = "atac", parallel = parallel, num.cores = cores)
  }
  
  #worker function
  atacCompartments <- function(obj, original.obj, res = 1e6, chr = NULL, targets = NULL,
                                genome = c("hg19", "hg38", "mm9", "mm10"),
                                prior.means = NULL, bootstrap = TRUE,
                                num.bootstraps = 1000, parallel = FALSE,
                                cores = 2, group = group, bootstrap.means = NULL) {
    #this is the main analysis function for computing compartments from atacs
    #make sure the input is sane
    if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")

    #what genome do we have
    genome <- match.arg(genome)
    
    #set the parallel back-end core number
    if (parallel) options(mc.cores = cores)
    
    #update
    message("Computing compartments for ", chr)
    obj <- keepSeqlevels(obj, chr, pruning.mode = "coarse")
    original.obj <- keepSeqlevels(original.obj, chr, pruning.mode = "coarse")
    
    #take care of the global means
    if (!is.null(prior.means)) {
      #this assumes that we've alread computed the global means
      pmeans <- as(prior.means, "GRanges")
      pmeans <- keepSeqlevels(pmeans, chr, pruning.mode = "coarse")
      #go back to a matrix
      prior.means <- as(pmeans, "matrix")
      colnames(prior.means) <- "globalMean"
    }
    
    #get the shrunken bins
    obj.bins <- shrinkBins(obj, original.obj, prior.means = prior.means, chr = chr,
                           res = res, targets = targets, assay = "atac",
                           genome = genome)
    #compute correlations
    if (group) obj.cor <- getCorMatrix(obj.bins, squeeze = FALSE)
    if (isFALSE(group)) obj.cor <- getCorMatrix(obj.bins, squeeze = TRUE)
    if (any(is.na(obj.cor$binmat.cor))) {
      obj.cor$gr$pc <- matrix(rep(NA, nrow(obj.cor$binmat.cor)))
      obj.svd <- obj.cor$gr
    } else {
      #compute SVD of correlation matrix
      obj.svd <- getABSignal(obj.cor, assay = "atac")
    }
    
    if (isFALSE(bootstrap)) return(obj.svd)
    
    #bootstrap the estimates
    #always compute confidence intervals too
    #take care of the global means
    if (bootstrap) {
      #this assumes that we've alread computed the global means
      bmeans <- as(bootstrap.means, "GRanges")
      bmeans <- keepSeqlevels(bmeans, chr, pruning.mode = "coarse")
      #go back to a matrix
      bmeans <- as(bmeans, "matrix")
      colnames(bmeans) <- rep("globalMean", ncol(bmeans))
    }
    
    obj.bootstrap <- bootstrapCompartments(obj, original.obj, bootstrap.samples = num.bootstraps,
                                           chr = chr, assay = "atac", parallel = parallel, cores = cores,
                                           targets = targets, res = res, genome = genome, q = 0.95,
                                           svd = obj.svd, group = group, bootstrap.means = bmeans)
    
    #combine and return
    return(obj.bootstrap)
  }
  
  #initialize global means
  #gmeans <- getGlobalMeans(obj, targets = targets, assay = "atac")
  
  if (parallel & isFALSE(group)) {
    atac.compartments <- mclapply(columns, function(s) {
      obj.sub <- obj[,s]
      message("Working on ", s)
      sort(unlist(as(lapply(chr, function(c) atacCompartments(obj.sub, obj, res = res,
                                                               chr = c, targets = targets, genome = genome,
                                                               bootstrap = bootstrap, prior.means = prior.means,
                                                               num.bootstraps = num.bootstraps, parallel = boot.parallel,
                                                               cores = boot.cores, group = group, bootstrap.means = bmeans)), "GRangesList")))
    }, mc.cores = cores)
  }
  
  if (!parallel & isFALSE(group)) {
    atac.compartments <- lapply(columns, function(s) {
      obj.sub <- obj[,s]
      message("Working on ", s)
      sort(unlist(as(lapply(chr, function(c) atacCompartments(obj.sub, obj, res = res,
                                                               chr = c, targets = targets, genome = genome,
                                                               bootstrap = bootstrap, prior.means = prior.means,
                                                               num.bootstraps = num.bootstraps, parallel = boot.parallel,
                                                               cores = boot.cores, group = group, bootstrap.means = bmeans)), "GRangesList")))
    })
  }
  
  if (parallel & isTRUE(group)) {
    atac.compartments <- sort(unlist(as(mclapply(chr, function(c) {
      atacCompartments(obj, obj, res = res,
                        chr = c, targets = targets, genome = genome,
                        bootstrap = bootstrap,num.bootstraps = num.bootstraps, prior.means = prior.means,
                        parallel = boot.parallel, cores = boot.cores, group = group, bootstrap.means = bmeans)}, mc.cores = cores),
      "GRangesList")))
  }
  
  if (!parallel & isTRUE(group)) {
    atac.compartments <- sort(unlist(as(lapply(chr, function(c) {
      atacCompartments(obj, obj, res = res,
                        chr = c, targets = targets, genome = genome, prior.means = prior.means,
                        bootstrap = bootstrap,num.bootstraps = num.bootstraps,
                        parallel = boot.parallel, cores = boot.cores, group = group, bootstrap.means = bmeans)}),
      "GRangesList")))
  }
  
  #if group-level treat a little differently
  if (group) {
    return(atac.compartments)
  }
  #convert to GRangesList
  atac.compartments <- as(atac.compartments, "CompressedGRangesList")
  #return as a RaggedExperiment
  return(RaggedExperiment(atac.compartments, colData = colData(obj)))
}
