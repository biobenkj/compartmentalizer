#' @title Fisher's Z transformation
#'
#' @description 
#' \code{fisherZ} returns (squeezed) Fisher's Z transformed Pearson's r
#'
#' @details 
#' This function returns (squeezed) Fisher's Z transformed Pearson's r
#'
#' @param cormat    Pearson correlation matrix
#'
#' @return    Fisher Z transformed Pearson correlations
#' @export
#'
#' @examples
#' 
#' #Generate a random binary (-1, 1) matrix
#' mat <- matrix(sample(c(1,-1), 10000, replace = TRUE), ncol = 100)
#' 
#' #Correct matrix diag
#' diag(mat) <- 1
#' 
#' #Transform
#' mat.transform <- fisherZ(mat)

fisherZ <- function(cormat) {
  if (any(cormat) == 1) cormat <- .squeezeit(cormat)
  atanh(cormat)
}

#' @title Fisher's Z transformation
#'
#' @description 
#' \code{fisherZ} returns the inverse (squeezed) Fisher's Z transformed Pearson's r. This will fail if a matrix is used as input instead of a vector.
#'
#' @details 
#' This function returns the inverse (squeezed) Fisher's Z transformed Pearson's r
#'
#' @param cormat    vector of Fisher's Z transformed Pearson correlations or an eignevector
#'
#' @return    Back transformed Fisher's Z
#' @export
#'
#' @examples
#' 
#' #Generate a random binary (-1, 1) matrix
#' mat <- matrix(sample(c(1,-1), 10000, replace = TRUE), ncol = 100)
#' 
#' #Correct matrix diag
#' diag(mat) <- 1
#' 
#' #Transform
#' mat.transform <- fisherZ(mat)
#' 
#' #Back transform
#' mat.transform.inverse <- apply(mat.transform, 1, ifisherZ)

ifisherZ <- function(cormat) {
  tanh(cormat)
}
