#' @title Estimate A/B compartments from WGBS data
#'
#' @description 
#' \code{getWGBSABsignal} returns estimated A/B compartments from whole genome bisulfite sequencing data
#'
#' @details 
#' This function estimates A/B compartments with some Dirichlet smoothing and shrinking towards a global mean of targets or across samples
#'
#' @param obj A BSseq object
#' @param res Compartment resolution (in bp)
#' @param globalMeanSet The global mean to shrink towards (calculated on the fly if NULL)
#' @param noMean Whether to shrink towards a global mean
#' @param targets Target samples to shrink towards
#' @param parallel Parallel estimation
#' @param allchrs Analyze all chromosomes (chr should be NULL if TRUE)
#' @param chr Chromosomes to analyze (list of arbitrary chromosomes)
#' @param regions GRanges object of loci locations (calculated on the fly)
#' @param genome Genome to analyze
#' @param preprocess Preprocessing with filtering, smoothing, and imputation
#' @param ... Other args
#'
#' @return compartment estimates
#' @import biscuiteer
#' @import bsseq
#' @import parallel
#' @import Homo.sapiens
#' @import impute
#' @export
#' 
#' @examples 
#' library(biscuiteer)
#' library(Homo.sapiens)
#' library(impute)
#' 
#' data(cell_cycle_hansen_chr22, package = "compartmap")
#' wgbs_compartments <- getWGBSABsignal(data.chr22, parallel = FALSE, chr = "chr22")

getWGBSABsignal <- function(obj, res=1e6, globalMeanSet = NULL, noMean = FALSE, targets = NULL, parallel=FALSE, allchrs=FALSE, chr = NULL, regions = NULL, genome = "hg19", preprocess = TRUE, ...) {
  
  # Run preprocessing if specified
  # Note: Dirichlet smoothing uses Jeffrey's prior as a value for k
  if (preprocess == TRUE) obj <- filterAndImputeWGBS(obj, k = 0.5, minCov = 3, minSamp = 2)
  
  columns <- colnames(obj)
  names(columns) <- columns
  
  if (any(is.na(obj))) stop("The input matrix cannot have any NA values in it. Either impute the missing values or remove loci/rows that contain NAs using filterAndImputeWGBS().")
  if (allchrs == TRUE & !is.null(chr)) stop("You've requested genome wide analysis in addition to specific chromosomes. Set chrs = 'all' for all chromosomes.")
  if (allchrs == FALSE & is.null(chr)) stop("You need to specify which chromosomes to use for analysis...")
  if (is.null(globalMeanSet)) globalMeanSet <- .getGlobalMeans(obj, targets)
  
  #Calculate things once... 
  if (is.null(regions)) regions <- .buildRanges(obj)
  
  #This assumes human... 
  #TODO: fix this to support more organisms...
  if (allchrs == TRUE) chr <- paste0("chr", c(seq(1,22)))
  
  #If we are squeezing to a global mean
  if (parallel & noMean == FALSE) {
    do.call(cbind, 
            mclapply(columns,.getWGBSShrunkenCompartments,obj=obj,globalMeanSet=globalMeanSet,noMean=noMean,targets=targets,regions=regions,genome=genome,chr=chr))
  } else if (noMean == FALSE) { 
    do.call(cbind, 
            lapply(columns,.getWGBSShrunkenCompartments,obj=obj,globalMeanSet=globalMeanSet,noMean=noMean,targets=targets,regions=regions,genome=genome,chr=chr))
  } else {
    #Not squeezing to global mean
    .getWGBSShrunkenCompartments(obj=obj,globalMeanSet=globalMeanSet,noMean=noMean,targets=targets,regions=regions,genome=genome,chr=chr)
  } 
}

#Helper function to find global means
.getGlobalMeans <- function(obj, targets = NULL) {
  if (!is.null(targets)){
    stargets <- .getShrinkageTargets(obj, targets)
    message(paste0("Using ", paste(shQuote(targets), collapse = ", "), " as shrinkage targets..."))
    meanBeta <- matrix(rowMeans(stargets, na.rm=TRUE), ncol=1)
  }
  else (meanBeta <- matrix(rowMeans(obj, na.rm=TRUE), ncol=1))
  colnames(meanBeta) <- "globalMean"
  return(meanBeta) 
}

#Helper function to identify a group of samples to use as shrinkage targets
#Subset a group of samples for binning
.getShrinkageTargets <- function(obj, group) {
  if (all(group %in% colnames(obj))) stargets.obj <- obj[,group]
  else (stop(paste0("Could not find ", group, " in the colnames of the input matrix...")))
  return(stargets.obj)
}


#A helper function to reconstruct the GRanges object for binning if not supplied from the input matrix
#These *should* be the rownames of the input matrix post biscuiteer
.buildRanges <- function(obj) {
  if (!is.null(rownames(obj))) {
    message("Constructing rowRanges from the input matrix...")
    snames <- sapply(strsplit(rownames(obj), ":"), `[`, 1)
    loc <- as.numeric(sapply(strsplit(rownames(obj), ":"), `[`, 2))
    myranges <- GRanges(seqnames = snames,
                        strand = "*",
                        ranges = IRanges(start = loc, width = 1))
  }
  else (stop("Could not reconstruct the rowRanges from the input matrix..."))
  return(myranges)
}

#Helper function to perform A/B compartment inference, shrinking towards a global mean (or not)
.getWGBSShrunkenCompartments <- function(column, obj, regions = NULL, globalMeanSet=NULL, noMean = FALSE, res=1e6, targets = NULL, chr = NULL, genome = "hg19", ...) {
  if (is.null(globalMeanSet)) globalMeanSet <- .getGlobalMeans(obj, targets)
  chrs <- chr
  names(chrs) <- chrs
  if (!is.null(globalMeanSet) & isFALSE(noMean)) {
    getMeanPairedChr <- function(chr) { 
      message("Computing shrunken eigenscores for ", column, " on ", chr, "...") 
      binmat <- getBinMatrix(cbind(obj[,column], globalMeanSet), regions, chr = chr, res = res, FUN = mean, genome = genome)
      cormat <- getCorMatrix(binmat, squeeze = TRUE)
      #Stupid check for perfect correlation with global mean
      if (any(is.na(cormat$binmat.cor))) {
        absig <- matrix(rep(NA, nrow(cormat$binmat.cor)))
      }
      else {
        absig <- getABSignal(cormat, squeeze = TRUE)$pc
      }
      
      return(absig)
    }
    unlist(lapply(chrs, getMeanPairedChr))
  }
  else if (isTRUE(noMean)) {
    getSmoothedPairedChr <- function(chr) { 
      message("Computing eigenscores on ", chr, "...") 
      binmat <- getBinMatrix(obj, regions, chr = chr, res = res, FUN = mean, genome = genome)
      cormat <- getCorMatrix(binmat)
      #Stupid check for NAs
      if (any(is.na(cormat$binmat.cor))) {
        warning("NAs were found in the correlation matrix. Eigenvector may not handle NAs.")
        absig <- try(getABSignal(cormat)$pc)
      }
      else {
        absig <- getABSignal(cormat)$pc
      }
      
      return(absig)
    }
    unlist(lapply(chrs, getSmoothedPairedChr))
  }
}
