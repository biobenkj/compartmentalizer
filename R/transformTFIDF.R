#' Transform/normalize compartment calls using TF-IDF
#' 
#' @name transformTFIDF
#'
#' @param obj n x p input matrix (n = samples/cells; p = compartments)
#' @param scale.factor Scaling factor for the term-frequency (TF)
#'
#' @return A TF-IDF transformed matrix of the same dimensions as the input
#'
#' @import Matrix
#'
#' @examples
#' 
#' m <- 1000
#' n <- 100
#' mat <- round(matrix(runif(m*n), m, n))
#' #Input needs to be a tall matrix
#' tfidf <- transformTFIDF(mat)
#' 
#' @export
transformTFIDF <- function(obj, scale.factor = 1e5) {
  #this filters using TF-IDF on a *matrix* object
  if (!is(obj, "matrix")) stop("Input needs to be a matrix.")
  #
  # FIXME: cite the following in the docs
  #
  #the following code was modeled after or taken from:
  #http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
  #https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/
  #binarize the matrix
  #this assumes n x p matrix (e.g. a wide matrix)
  #check and transpose as needed
  #input matrix is tall
  if (dim(obj)[1] > dim(obj)[2]) obj <- t(obj)
  #make sparse
  obj.binary <- Matrix(.binarizeMatrix(t(obj)), sparse = TRUE)
  #compute term-frequency
  tf <- t(t(obj.binary) / Matrix::colSums(obj.binary))
  #scale
  tf@x <- log1p(tf@x * scale.factor)
  #inverse-document frequency smooth
  idf <- log(1 + ncol(obj.binary) / Matrix::rowSums(obj.binary))
  #transform
  tfidf <- .tfidf(tf, idf)
  #cast back to a matrix since things like UMAP don't like sparse matrices
  tfidf <- as.matrix(tfidf)
  return(t(tfidf))
}


#helper function
.binarizeMatrix <- function(obj) {
  #set positive values to 1 and negative to 0
  #open chromatin in atac or RNA is 1
  #closed chromatin in bisulfite or arrays is 1
  #just associate 1 with signal
  obj[obj > 0] <- 1
  obj[obj < 0] <- 0
  return(obj)
}

#helper function for TF-IDF transform
#modeled after http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
.tfidf <- function(tf, idf) {
  tf = t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = t(tf)
  return(tf)
}

#' Transform/normalize compartment calls using TF-IDF on HDF5-backed objects
#' 
#' @name hdf5TFIDF
#'
#' @param h5 SummarizedExperiment object, DelayedMatrix, or a normal matrix
#' @param scale.factor Scaling factor for the term-frequency (TF)
#' @param return.dense Whether to return a dense, in memory matrix
#' @param return.se Whether to return the TF-IDF matrix as a new assay in the SummarizedExperiment
#'
#' @return A TF-IDF transformed matrix of the same dimensions as the input
#' @import Matrix
#' @import DelayedMatrixStats
#' @import DelayedArray
#' @import HDF5Array
#'
#' @examples
#' 
#' m <- 1000
#' n <- 100
#' mat <- round(matrix(runif(m*n), m, n))
#' #Input needs to be a tall matrix
#' tfidf <- hdf5TFIDF(mat)
#' 
#' @export
hdf5TFIDF <- function(h5, scale.factor = 1e5,
                      return.dense = FALSE,
                      return.se = FALSE) {
  #binarze
  if (is(h5, "SummarizedExperiment")) {
    if (!is(assay(h5), "DelayedMatrix")) {
      #coerce to hdf5 backing to work with any SE
      assay(h5) <- as(assay(h5), "HDF5Matrix")
    }
    assay(h5)[assay(h5) > 0] <- 1
    #tall matrix
    h5.mat <- assay(h5)
  }
  if (is(h5, "DelayedMatrix")) {
    #make the matrix tall if needed
    if (dim(h5)[1] < dim(h5)[2]) h5 <- t(h5)
    h5[h5 > 0] <- 1
    h5.mat <- h5
  }
  if (is(h5, "matrix")) {
    h5[h5 > 0] <- 1
    h5.mat <- as(h5, "HDF5Matrix")
  }
  #term frequency
  message("Computing term frequency.")
  tf <- t(t(h5.mat)/DelayedMatrixStats::colSums2(h5.mat))
  #scale
  tf <- log1p(tf * scale.factor)
  #inverse document frequency
  message("Computing inverse document frequency.")
  idf <- log(1 + ncol(h5.mat)/DelayedMatrixStats::rowSums2(h5.mat))
  #cast the tf matrix back to a sparse matrix
  #TODO: fix this ugliness...
  tf.mat <- as.matrix(tf)
  tf.sparse <- Matrix(tf.mat, sparse = TRUE)
  #transpose for TF-IDF
  tf.sparse <- t(tf.sparse)
  #TF-IDF applied
  message("TF-IDF")
  tf.sparse@x <- tf.sparse@x * rep.int(idf, diff(tf.sparse@p))
  #transpose again
  tf.sparse <- t(tf.sparse)
  #coerce back to dense matrix
  if (return.dense) {
    message("WARNING: This might blow up!")
    message("If you get a cholmod error: problem too large, set return.dense to FALSE.")
    message("You will get a sparse matrix returned instead.")
    return(as.matrix(tf.sparse))
  }
  if (return.se) {
    if (!is(h5, "SummarizedExperiment")) {
      return(tf.sparse)
    }
    message("Returning the TF-IDF matrix into the SummarizedExperiment.")
    assays(h5)$tfidf <- tf.sparse
    return(h5)
  }
  return(tf.sparse)
}
