#' Compute the sparsity of each bin that will be used to infer compartments
#'
#' @name getBinSparsity
#' 
#' @param mat Input matrix (samples as columns and loci as rows)
#' @param row.ranges The rowRanges or a GRanges object corresponding to genomic locations of loci
#' @param chr Which chromosome to operate on
#' @param res What resolution or bin size to use (default: 1e5)
#' @param min.depth Minimum depth to consider (default: 1)
#' @param scale.by.nonzero.bins Whether to scale the sparsity measure by the number of non-zero bins (default: TRUE)
#' @param filter.zero.bins Whether to filter non-zero bins (default: TRUE)
#' @param scale.by.resolution Whether to scale the per-locus sparsity by the resolution (default: FALSE)
#' @param species What species to use to determine chromosome sizes
#' 
#' @import GenomicRanges
#' @import pbapply
#' @import Homo.sapiens
#' @import Mus.musculus
#' @import SummarizedExperiment
#'
#' @return A list object containing either the scaled/summarized sparsity measures or a per-locus sparsity matrix and a corresponding GRanges object
#' @export
#'
#' @examples
#'

getBinSparsity <- function(mat, row.ranges, chr = "chr14",
                           res = 1e5, min.depth = 1, scale.by.nonzero.bins = TRUE,
                           filter.zero.bins = TRUE, scale.by.resolution = FALSE,
                           species = c("human", "mouse")) {
  #get signal sparsity within each bin
  #this could be used to filter or plot signal sparsity in samples
  #use this to filter out cells of poor quality by scaling by number of bins
  #e.g. colSums(sparsity.mat) / bins | resolution across the genome
  #this makes sense because it's signal per bin
  #expects input to be a matrix, much like getBinMatrix()
  #most of this code is duplicated from getBinMatrix()
  
  if (any(is.na(mat))){
    stop("Matrix must not contain NAs")
  }
  if (nrow(mat)!=length(row.ranges)){
    stop("Provided GRanges must have length equal to the matrix number of rows")
  }
  
  #figure out which species we are using
  species <- match.arg(species)
  message("Using ", species, " as the target species.")
  
  message("Working on ", chr)
  if (species == "human") {
    chr.end <- seqlengths(Homo.sapiens)[chr]
  } else {
      chr.end <- seqlengths(Mus.musculus)
    }
  chr.start <- 0 #0-based indexing
  
  start <- seq(chr.start, chr.end, res) #Build the possible bin ranges given resolution
  end <- c(start[-1], chr.end) - 1L #Set the end ranges for desired resolution
  
  #Build up the genomic ranges object given chr, start, end, and resolution
  gr.bin <- GRanges(seqnames = chr,
                    ranges = IRanges(start = start, end = end))
  
  #Identify overlaps between the user defined GRanges object (loci) and bins
  ids <- findOverlaps(row.ranges, gr.bin, select="first")
  
  #Get the number of bins overlapping loci
  n <- length(gr.bin)
  message(n, " bins created.")
  
  signalSparsity <- function(bin, min.depth = min.depth) {
    ss <- table(bin >= min.depth)["TRUE"]
    return(ifelse(is.na(ss), 0, ss))
  }
  
  #summarize signal density within each bin, given resolution
  message("Computing sparsity measures.")
  sparsity.mat <- pbapply(mat, 2, function(f) {
    zvec <- rep(0, n) #Generate a vector of zeroes
    a <- tapply(f, INDEX=ids, FUN=signalSparsity, min.depth = min.depth) #Summarize data
    zvec[as.numeric(names(a))] <- a
    return(zvec)
  })
  
  colnames(sparsity.mat) <- colnames(mat) #Set colnames
  
  #filter out bins without any signal in them
  if (filter.zero.bins) {
    message("Filtering out zero count bins across samples.")
    nonzero.bins <- rowSums(sparsity.mat) != 0
    sparsity.mat <- sparsity.mat[nonzero.bins,]
    gr.bin <- gr.bin[nonzero.bins,]
    message("Filtered ", table(nonzero.bins)["FALSE"], " bins from the data.")
  }
  
  #scale the per-bin signal by the resolution
  if (scale.by.resolution) {
    message("Scaling the per-bin signal by the resolution:", res, ".")
    sparsity.mat <- log((sparsity.mat / res) + 1)
  }
  
  #return scaled sparsity measure
  if (scale.by.nonzero.bins) {
    message("Scaling sparsity measure.")
    sparsity.mat <- colSums(sparsity.mat) / nrow(sparsity.mat)
    return(list(scaled.sparsity.mat=sparsity.mat,
                binned.granges=gr.bin))
  }

  return(list(sparsity.mat=sparsity.mat,
              binned.granges=gr.bin))
}
