#' @title Estimate A/B compartments from BS-seq data
#'
#' @description 
#' \code{getBSABsignal} returns estimated A/B compartments from BS-seq data.
#'
#' @param obj Input SummarizedExperiment object
#' @param res Compartment resolution in bp
#' @param parallel Whether to run samples in parallel
#' @param chr What chromosome to work on (leave as NULL to run on all chromosomes)
#' @param targets Samples/cells to shrink towards
#' @param preprocess Whether to preprocess the bs-seq data prior to compartment inference
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
#' data("meth_wgbs_450k_chr14", package = "compartmap")
#' wgbs_compartments <- getWGBSABsignal(wgbs.data.chr14, parallel=F, chr="chr14", bootstrap=F, genome="hg19", wgbs.type="hm450")

getBSABsignal <- function(obj, res = 1e6, parallel = TRUE, chr = NULL,
                             targets = NULL, preprocess = TRUE, cores = 2,
                             bootstrap = TRUE, num.bootstraps = 1000,
                             genome = c("hg19", "hg38", "mm9", "mm10"),
                             other = NULL,
                             group = FALSE, boot.parallel = TRUE, boot.cores = 2) {
  
  #preprocess the wgbss
  if (preprocess) {
    obj <- preprocessBSseq(obj = obj, genome = genome, other = other)
  }
  
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
  
  #worker function
  wgbsCompartments <- function(obj, original.obj, res = 1e6, chr = NULL, targets = NULL,
                                genome = c("hg19", "hg38", "mm9", "mm10"),
                                prior.means = NULL, bootstrap = TRUE,
                                num.bootstraps = 1000, parallel = FALSE,
                                cores = 2, group = group) {
    #this is the main analysis function for computing compartments from wgbss
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
    
    #get the shrunken bins
    obj.bins <- shrinkBins(obj, original.obj, prior.means = prior.means, chr = chr,
                           res = res, targets = targets, assay = "bisulfite",
                           genome = genome)
    #compute correlations
    obj.cor <- getCorMatrix(obj.bins, squeeze = TRUE)
    if (any(is.na(obj.cor$binmat.cor))) {
      obj.cor$gr$pc <- matrix(rep(NA, nrow(obj.cor$binmat.cor)))
      obj.svd <- obj.cor$gr
    } else {
      #compute SVD of correlation matrix
      obj.svd <- getABSignal(obj.cor, assay = "bisulfite")
    }
    
    if (isFALSE(bootstrap)) return(obj.svd)
    
    #bootstrap the estimates
    #always compute confidence intervals too
    obj.bootstrap <- bootstrapCompartments(obj, original.obj, bootstrap.samples = num.bootstraps,
                                           chr = chr, assay = "bisulfite", parallel = parallel, cores = cores,
                                           targets = targets, res = res, genome = genome, q = 0.95,
                                           svd = obj.svd, group = group)
    
    #combine and return
    return(obj.bootstrap)
  }
  
  #initialize global means
  #gmeans <- getGlobalMeans(obj, targets = targets, assay = "bisulfite")
  
  if (parallel & isFALSE(group)) {
    wgbs.compartments <- mclapply(columns, function(s) {
      obj.sub <- obj[,s]
      message("Working on ", s)
      sort(unlist(as(lapply(chr, function(c) wgbsCompartments(obj.sub, obj, res = res,
                                                               chr = c, targets = targets, genome = genome,
                                                               bootstrap = bootstrap,
                                                               num.bootstraps = num.bootstraps, parallel = boot.parallel,
                                                               cores = boot.cores, group = group)), "GRangesList")))
    }, mc.cores = cores)
  }
  
  if (!parallel & isFALSE(group)) {
    wgbs.compartments <- lapply(columns, function(s) {
      obj.sub <- obj[,s]
      message("Working on ", s)
      sort(unlist(as(lapply(chr, function(c) wgbsCompartments(obj.sub, obj, res = res,
                                                               chr = c, targets = targets, genome = genome,
                                                               bootstrap = bootstrap,
                                                               num.bootstraps = num.bootstraps, parallel = boot.parallel,
                                                               cores = boot.cores, group = group)), "GRangesList")))
    })
  }
  
  if (parallel & isTRUE(group)) {
    wgbs.compartments <- sort(unlist(as(mclapply(chr, function(c) {
      wgbsCompartments(obj, obj, res = res,
                        chr = c, targets = targets, genome = genome,
                        bootstrap = bootstrap,num.bootstraps = num.bootstraps,
                        parallel = boot.parallel, cores = boot.cores, group = group)}, mc.cores = cores),
      "GRangesList")))
  }
  
  if (!parallel & isTRUE(group)) {
    wgbs.compartments <- sort(unlist(as(lapply(chr, function(c) {
      wgbsCompartments(obj, obj, res = res,
                        chr = c, targets = targets, genome = genome,
                        bootstrap = bootstrap,num.bootstraps = num.bootstraps,
                        parallel = boot.parallel, cores = boot.cores, group = group)}),
      "GRangesList")))
  }
  
  #if group-level treat a little differently
  if (group) {
    return(wgbs.compartments)
  }
  #convert to GRangesList
  wgbs.compartments <- as(wgbs.compartments, "CompressedGRangesList")
  #return as a RaggedExperiment
  return(RaggedExperiment(wgbs.compartments, colData = colData(obj)))
}

#' Preprocess BS-seq for compartment inference
#'
#' @name preprocessBSseq
#'
#' @param obj Input SummarizedExperiment
#' @param genome What genome are we working on ("hg19", "hg38", "mm9", "mm10")
#' @param other Another arbitrary genome to compute compartments on
#'
#' @return A preprocessed SummarizedExperiment to compute compartments
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 
preprocessBSseq <- function(obj,
                             genome = c("hg19", "hg38", "mm9", "mm10"),
                             other = NULL) {
  #make sure the input is sane
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment")
  
  #what genome do we have
  genome <- match.arg(genome)
  
  #smooth
  obj.smooth <- smoothBSCounts(obj)
  
  #subset the wgbs to open sea CpGs
  obj.opensea <- filterOpenSea(obj.smooth, genome = genome, other = other)
  
  #convert things to M-values
  #check the names of the assays
  if (!any(getAssayNames(obj.opensea) %in% c("counts"))) {
    stop("The assays slot should contain 'counts' for bisulfite data.")
  }
  
  #impute missing values and drop samples that are too sparse
  message("Imputing missing values.")
  obj.opensea.imputed <- imputeKNN(obj.opensea, assay = "bisulfite")
  
  return(obj.opensea.imputed)
}
