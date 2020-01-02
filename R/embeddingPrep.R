#' Filter and convert/cluster compartment calls for embedding and dimension reduction
#'
#' @name embeddingPrep
#'
#' @param obj Input either a RaggedExperiment or output from fixCompartments
#' @param filter Whether to filter compartments
#' @param min.conf Minimum confidence estimate to use when filtering
#' @param min.eigen Minimum absolute eigenvalue to use when filtering
#' @param tf.idf Whether to TF-IDF transform the data
#' @param cluster Whether to identify clusters in the data
#' @param k How many neighbors when identifying clusters
#' @param dims How many dimensions to use for identifying clusters
#'
#' @return A prepared matrix of list of compartment calls for embedding and dimension reduction
#' @import RaggedExperiment
#' @import SummarizedExperiment
#' @import Matrix
#' @import BiocSingular
#' @import RANN
#' @import igraph
#' @export
#'
#' @examples
#' 
embeddingPrep <- function(obj, filter = TRUE,
                          min.conf = 0.7, min.eigen = 0.02,
                          tf.idf = FALSE, cluster = TRUE,
                          k = 10, dims = NULL) {
  
  #check and see if we have just one sample and stop if we do
  if (is(obj, "GRanges")) stop("Need at least 2 samples.")
  
  #fix up the compartments
  if (is(obj, "RaggedExperiment")) obj <- fixCompartments(obj, min.conf = min.conf)

  #filter compartments down
  if (isTRUE(filter)) {
    obj <- filterCompartments(obj,
                              min.conf = min.conf,
                              min.eigen = min.eigen)
  }
  
  #set a seed
  set.seed(1000)
  #turn into a matrix for embedding approaches
  #do a little RaggedExperiment magic
  #to allow for disjoint compartment calls
  if (is(obj, "list")) obj <- RaggedExperiment(as(obj, "GRangesList"))
  
  #condense
  fixed <- ifelse("flip.score" %in% getAssayNames(obj),
                  TRUE, FALSE)
  if (isTRUE(fixed)) se <- compactSummarizedExperiment(obj, i = "flip.score")
  if (isFALSE(fixed)) se <- compactSummarizedExperiment(obj, i = "score")
  
  #go to a matrix
  score.mat <- t(as.matrix(assay(se)))
  
  #make NAs soft zeros
  message("Making NAs soft zeros.")
  score.mat[is.na(score.mat)] <- 0
  
  #tf-idf transform
  if (isTRUE(tf.idf)) {
    message("Transforming with TF-IDF.")
    score.mat <- transformTFIDF(score.mat)
  }
  
  if (isFALSE(cluster)) return(score.mat)
  
  #dim reduction
  if (!is.null(dims) && dims <= nrow(score.mat)) {
    message("Running PCA and retaining ", dims, " dimensions.")
    pca.results <- runPCA(t(score.mat), rank = dims, BSPARAM=ExactParam())
    rownames(pca.results$rotation) <- rownames(score.mat)
  } 
  if (!is.null(dims) && dims > nrow(score.mat)) {
    stop("dims needs to be less than the number of samples/cells we have.")
    }
  
  #cluster
  message("Clustering using k=", k, " neighbors.")
  if (exists("pca.results")) {
    #compute neighbors in PCA space
    knn <- RANN::nn2(pca.results$rotation, k = k)
    adj_mat <- .getAdjacencyMatrix(pca.results$rotation, knn)
  } else {
    #compute neighbors in TF-IDF or un-transformed space
    knn <- RANN::nn2(score.mat, k = k)
    adj_mat <- .getAdjacencyMatrix(score.mat, knn)
  }
  
  #convert to graph
  adjacency_graph <- igraph::graph.adjacency(adj_mat,
                                             mode = "undirected")
  adjacency_graph <- igraph::simplify(adjacency_graph)
  
  #find communities/clusters
  message("Finding communities in the graph.")
  louvain_com <- igraph::cluster_louvain(adjacency_graph)
  com <- louvain_com$membership
  names(com) <- louvain_com$names
  if (isTRUE(cluster) && !is.null(dims)) return(list(compartment.mat=score.mat,
                                                     communities=com,
                                                     reduced_dims=pca.results))
  if (isTRUE(cluster)) return(list(compartment.mat=score.mat,
                                   communities=com))
}

#' Filter compartments using confidence estimates and eigenvalue thresholds
#' 
#' @name filterCompartments
#'
#' @param obj Output of condenseSE or fixCompartments
#' @param min.conf Minimum confidence estimate to use when filtering
#' @param min.eigen Minimum absolute eigenvalue to use when filtering
#'
#' @return A filtered/subset of the input object/list
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 
filterCompartments <- function(obj, min.conf = 0.7, min.eigen = 0.02) {
  #filter compartments
  message("Filtering compartments based on a minimum confidence of ", min.conf*100, "%")
  message("Filtering compartments based on a minimum absolute eigen value of ", min.eigen)
  if (is(obj, "list")) {
    filt.compartments <- lapply(obj, function(x) {
      filt <- apply(mcols(x), 1, function(r) {
        #check if we have "fixed" things
        if ("flip.score" %in% names(r)) {
          filt.score <- ifelse(as.numeric(r["flip.conf.est"]) >= min.conf &
                                 abs(as.numeric(r["flip.score"])) >= min.eigen,
                               TRUE, FALSE)
          return(filt.score)
        } else {
          return(ifelse(as.numeric(r["conf.est"]) >= min.conf &
                          abs(as.numeric(r["score"])) >= min.eigen, TRUE, FALSE))
        }
      })
      return(x[filt,])
    })
  } else {
    filt <- apply(mcols(obj), 1, function(r) {
      #check if we have "fixed" things
      if ("flip.score" %in% names(r)) {
        filt.score <- ifelse(as.numeric(r["flip.conf.est"]) >= min.conf &
                               abs(as.numeric(r["flip.score"])) >= min.eigen,
                             TRUE, FALSE)
        return(filt.score)
      } else {
        return(ifelse(as.numeric(r["conf.est"]) >= min.conf &
                        abs(as.numeric(r["score"])) >= min.eigen, TRUE, FALSE))
      }
    })
    return(obj[filt,])
  }
}

#' Transform/normalize compartment calls using TF-IDF
#' 
#' @name transformTFIDF
#'
#' @param obj n x p input matrix (n = samples/cells; p = compartments)
#' @param scale.factor Scaling factor for the term-frequency (TF)
#'
#' @return A TF-IDF transformed matrix of the same dimensions as the input
#' @import Matrix
#' @export
#'
#' @examples
#' 
transformTFIDF <- function(obj, scale.factor = 1e5) {
  #this filters using TF-IDF on a *matrix* object
  if (!is(obj, "matrix")) stop("Input needs to be a matrix.")
  #the following code was modeled after or taken from:
  #http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
  #https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/
  #binarize the matrix
  #this assumes n x p matrix (e.g. a wide matrix)
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

#helper function for RANN results to adjacency matrix
#modeled from #https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/
.getAdjacencyMatrix <- function(obj, knn) {
  #knn input is KNN results
  knn_map <- knn$nn.idx
  #create an empty matrix
  adjacency_mat <- matrix(0, nrow(obj), nrow(obj))
  #fill it in
  rownames(adjacency_mat) <- colnames(adjacency_mat) <- rownames(obj)
  for(i in seq_len(nrow(obj))) {
    adjacency_mat[i,rownames(obj)[knn_map[i,]]] <- 1
  }
  return(adjacency_mat)
}
