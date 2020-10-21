#' Filter to open sea CpG loci
#'
#' @name filterOpenSea
#'
#' @param obj Input SummarizedExperiment or GRanges object
#' @param genome Which genome to filter
#'
#' @return Filtered to open sea CpG loci
#' @import SummarizedExperiment
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' opensea <- filterOpenSea(array.data.chr14, genome = "hg19")
#' 
#' @export
filterOpenSea <- function(obj, genome = c("hg19", "hg38", "mm10", "mm9"), other = NULL) {
  #get the desired open sea loci given the genome
  genome <- match.arg(genome)
  if (is.null(other)) {
    openseas.genome <- switch(genome,
                              hg19=data("openSeas.hg19", package="compartmap"),
                              hg38=data("openSeas.hg38", package="compartmap"),
                              mm10=data("openSeas.mm10", package="compartmap"),
                              mm9=data("openSeas.mm9", package="compartmap"))
  } else {
    #check if it's a GRanges flavored object
    if (!is(other, "GRanges")) stop("The 'other' input needs to be a GRanges of open sea regions")
    openseas.genome <- other
  }
  #Subset by overlaps
  message("Filtering to open sea CpG loci...")
  obj.openseas <- subsetByOverlaps(obj, get(openseas.genome))
  return(obj.openseas)
}

#' Gather open sea CpG from a GRanges of CpG islands
#' 
#' @description This function accepts a GRanges input of CpG islands that can
#' be derived from UCSC table browser and rtracklayer::import(yourbed.bed)
#'
#' @name filterOpenSea
#'
#' @param gr Input GRanges of CpG islands
#'
#' @return GRanges object that can be used with filterOpenSea()
#' @import rtracklayer
#' @import GenomicRanges
#' @export
#'
#' @examples
#' cpgi <- rtracklayer::import(system.file("inst/extdata/mm10_cpgi.bed", package = "compartmap"))
#' opensea_cpg <- getOpenSeas(cpgi)
#' 

getOpenSeas <- function(gr) {
  resorts <- trim(resize(gr, width(gr) + 8000, fix = "center"))
  openSeas <- subset(gaps(resorts), strand == "*")
  return(openSeas)
}

#' Get the open and closed compartment calls based on sign of singular values
#'
#' @param gr Input GRanges with associated mcols that represent singular values
#' @param cutoff Threshold to define open and closed states
#' @param assay The type of assay we are working with
#'
#' @return A vector of binary/categorical compartment states
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 
#' dummy <- matrix(rnorm(10000), ncol=25)
#' set.seed(1000)
#' sing_vec <- getSVD(dummy, k = 1, sing.vec = "left")
#' 
extractOpenClosed <- function(gr, cutoff = 0,
                              assay = c("array", "atac", "bisulfite")){
  #check for input to be GRanges
  if (!is(gr, "GRanges")) stop("Input needs to be a GRanges.")
  if (!("pc" %in% names(mcols(gr)))) stop("Need to have an mcols column be named 'pc'.")
  assay <- match.arg(assay)
  if (assay %in% c("array", "bisulfite")) return(ifelse(gr$pc < cutoff, "open", "closed"))
  if (assay == "atac") return(ifelse(gr$pc < cutoff, "closed", "open"))
}

#' Check if the assay is a SummarizedExperiment
#'
#' @param obj Input object
#'
#' @return Boolean
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' checkAssayType(array.data.chr14)

checkAssayType <- function(obj) {
  #helper function to check the class of an object
  is(obj, "SummarizedExperiment")
}

#' Get the assay names from a SummarizedExperiment object
#'
#' @param se Input SummarizedExperiment object
#'
#' @return The names of the assays
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' getAssayNames(array.data.chr14)

getAssayNames <- function(se) {
  #helper function to check the assay slot names
  names(assays(se))
}

#' Remove rows with NAs exceeding a threshold
#'
#' @param se Input SummarizedExperiment object
#' @param rowmax The maximum NAs allowed in a row as a fraction
#' @param assay The type of assay we are working with
#'
#' @return A filtered matrix
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' cleanAssayRows(array.data.chr14, assay = "array")

cleanAssayRows <- function(se, rowmax = 0.5,
                           assay = c("array", "atac", "bisulfite")) {
  assay <- match.arg(assay)
  switch(assay,
         array = se[rowMeans(is.na(assays(se)$Beta)) < rowmax,],
         atac = se[rowMeans(is.na(assays(se)$counts)) < rowmax,],
         bisulfite = se[rowMeans(is.na(assays(se)$counts)) < rowmax,])
}

#' Remove columns/cells/samples with NAs exceeding a threshold
#'
#' @param se Input SummarizedExperiment object
#' @param colmax The maximum number of NAs allowed as a fraction
#' @param assay The type of assay we are working with
#'
#' @return A filtered matrix
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' cleanAssayCols(array.data.chr.14, assay = "array")
cleanAssayCols <- function(se, colmax = 0.8,
                           assay = c("array", "atac", "bisulfite")) {
  assay <- match.arg(assay)
  switch(assay,
         array = se[,colMeans(is.na(assays(se)$Beta)) < colmax],
         atac = se[,colMeans(is.na(assays(se)$counts)) < colmax],
         bisulfite = se[,colMeans(is.na(assays(se)$counts)) < colmax])
}

#' Helper function: squeezed logit
#'
#' @param p       a vector of values between 0 and 1 inclusive
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between -Inf and +Inf
#'
#' @examples
#'
#'   set.seed(1234)
#'   p <- runif(n=1000)
#'   summary(p) 
#' 
#'   sqz <- 1 / (10**6)
#'   x <- flogit(p, sqz=sqz)
#'   summary(x) 
#'
#'   all( abs(p - fexpit(x, sqz=sqz)) < sqz )
#'   all( abs(p - fexpit(flogit(p, sqz=sqz), sqz=sqz)) < sqz ) 
#'
#' @export 
flogit <- function(p, sqz=0.000001) { 

  midpt <- 0.5
  deflate <- 1 - (sqz * midpt)
  if (any(p > 1 | p < 0, na.rm = TRUE)) stop("Values of p outside (0,1) detected.")
  squoze <- ((p - midpt) * deflate) + midpt
  return( log( squoze / (1 - squoze)) )
  
}

#' Helper function: expanded expit
#'
#' @param x       a vector of values between -Inf and +Inf
#' @param sqz     the amount by which we 'squoze', default is .000001
#'
#' @return        a vector of values between 0 and 1 inclusive
#'
#' @examples
#'
#'   set.seed(1234)
#'   x <- rnorm(n=1000)
#'   summary(x) 
#'
#'   sqz <- 1 / (10**6)
#'   p <- fexpit(x, sqz=sqz)
#'   summary(p)
#'
#'   all( (abs(x - flogit(p)) / x) < sqz )
#'   all( abs(x - flogit(fexpit(x))) < sqz )
#'
#' @export 
fexpit <- function(x, sqz=0.000001) {

  midpt <- .5
  squoze <- exp(x)/(1 + exp(x))
  inflate <- 1 / (1 - (sqz * midpt))
  p <- ((squoze - midpt) * inflate) + midpt 
  return(p)

}

#' Get the chromosomes from an object
#'
#' @param obj Input SummarizedExperiment object
#'
#' @return A character vector of chromosomes present in an object
#' @import SummarizedExperiment
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' getChrs(array.data.chr14)
#' 
#' @export
getChrs <- function(obj) {
  #get the chromosomes present in the object
  return(unique(as.character(seqnames(obj))))
}

#' Remove bootstrap estimates that failed
#'
#' @param obj Input list object with elements 'pc' and 'gr'
#'
#' @return A filtered list object
#' @export
removeEmptyBoots <- function(obj) {
  #remove NAs from a bootstrap list
  #this can happen if the correlation between the bins and eigenvector fails
  #theoretically we can recover these but need an additional utility to find consensus
  na.filt <- unlist(lapply(obj, function(n) ifelse(any(is.na(n)), FALSE, TRUE)))
  obj <- obj[na.filt]
  return(obj)
}

#' Get the seqlengths of a chromosome
#' 
#' The goal for this function is to eliminate the need to lug around
#' large packages when we only want seqlengths for things.
#'
#' @param genome The desired genome to use ("hg19", "hg38", "mm9", "mm10")
#' @param chr What chromosome to extract the seqlengths of
#'
#' @return The seqlengths of a specific chromosome
#' @import GenomicRanges
#'
#' @examples
#' hg19.chr14.seqlengths <- getSeqLengths(genome = "hg19", chr = "chr14")
#'
#' @export
getSeqLengths <- function(genome = c("hg19", "hg38", "mm9", "mm10"),
                          chr = "chr14") {
  #eventually we should support arbitrary genomes
  genome <- match.arg(genome)
  #check if the genome used exists in what is currently supported, stopping if not
  if (!genome %in% c("hg19", "hg38", "mm9", "mm10")) stop("Only human and mouse are supported for the time being.")
  #import
  genome.gr <- switch(genome,
                      hg19 = data("hg19.gr", package = "compartmap"),
                      hg38 = data("hg38.gr", package = "compartmap"),
                      mm9 = data("mm9.gr", package = "compartmap"),
                      mm10 = data("mm10.gr", package = "compartmap"))
  #make sure that the chromosome specified exists in the seqlevels
  if (!chr %in% seqlevels(get(genome.gr))) stop("Desired chromosome is not found in the seqlevels of ", genome)
  #get the seqlengths
  sl <- seqlengths(get(genome.gr))[chr]
  return(sl)
}
  
#' Get chunked sets of row-wise or column-wise indices of a matrix 
#' 
#' @name getMatrixBlocks
#'
#' @param mat Input matrix
#' @param chunk.size The size of the chunks to use for coercion
#' @param by.row Whether to chunk in a row-wise fashion
#' @param by.col Whether to chunk in a column-wise fashion
#'
#' @return A set of chunked indices
#'
#' @examples
#' #make a sparse binary matrix
#' library(Matrix)
#' m <- 100
#' n <- 1000
#' mat <- round(matrix(runif(m*n), m, n))
#' mat.sparse <- Matrix(mat, sparse = TRUE)
#' 
#' #get row-wise chunks of 10
#' chunks <- getMatrixBlocks(mat.sparse, chunk.size = 10)
#'
#' @export
getMatrixBlocks <- function(mat, chunk.size = 1e5,
                            by.row = TRUE, by.col = FALSE) {
  message("Using chunk size: ", chunk.size)
  if (by.row) {
    message("Breaking into row chunks.")
    return(split(1:nrow(mat), ceiling(seq_along(1:nrow(mat))/chunk.size)))
  }
  
  #assumes column-wise chunking
  message("Breaking into column chunks.")
  return(split(1:ncol(mat), ceiling(seq_along(1:ncol(mat))/chunk.size)))
}

#' Convert a sparse matrix to a dense matrix in a block-wise fashion 
#' 
#' @name sparseToDenseMatrix
#'
#' @param mat Input sparse matrix
#' @param blockwise Whether to do the coercion in a block-wise manner
#' @param by.row Whether to chunk in a row-wise fashion
#' @param by.col Whether to chunk in a column-wise fashion
#' @param chunk.size The size of the chunks to use for coercion
#' @param parallel Whether to perform the coercion in parallel
#' @param cores The number of cores to use in the parallel coercion
#'
#' @return A dense matrix of the same dimensions as the input
#' 
#' @import Matrix
#' @import parallel
#' 
#' 
#' @examples
#' #make a sparse binary matrix
#' library(Matrix)
#' m <- 100
#' n <- 1000
#' mat <- round(matrix(runif(m*n), m, n))
#' mat.sparse <- Matrix(mat, sparse = TRUE)
#' 
#' #coerce back
#' mat.dense <- sparseToDenseMatrix(mat.sparse, chunk.size = 10)
#' 
#' #make sure they are the same dimensions
#' dim(mat) == dim(mat.dense)
#' 
#' #make sure they are the same numerically
#' all(mat == mat.dense)
#'
#' @export 
sparseToDenseMatrix <- function(mat, blockwise = TRUE,
                                by.row = TRUE, by.col = FALSE,
                                chunk.size = 1e5, parallel = FALSE,
                                cores = 2) {
  if (isFALSE(blockwise)) return(as(mat, "matrix"))
  
  #do block-wise reconstruction of matrix
  chunks <- getMatrixBlocks(mat, chunk.size = chunk.size,
                            by.row = by.row, by.col = by.col)
  
  if (by.row & parallel) {
    return(do.call("rbind", mclapply(chunks, function(r) {
      return(as(mat[r,], "matrix"))
    }, mc.cores = cores)))
  }
  
  if (by.row & !parallel) {
    return(do.call("rbind", lapply(chunks, function(r) {
      return(as(mat[r,], "matrix"))
    })))
  }
  
  #assumes column-wise conversion
  if (by.col & parallel) {
    return(do.call("cbind", mclapply(chunks, function(r) {
      return(as(mat[,r], "matrix"))
    }, mc.cores = cores)))
  }
  
  return(do.call("cbind", lapply(chunks, function(r) {
    return(as(mat[,r], "matrix"))
  })))
  
}  
