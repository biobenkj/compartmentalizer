#' Compute the SVD of a matrix using irlba
#'
#' @name getSVD
#'
#' @param matrix A p x n input matrix
#' @param k Number of singular vectors to return
#'
#' @return A singular vector or matrix with sign corresponding to positive values
#' 
#' @import BiocSingular
#' @export
#'
#' @examples
#' 
#' dummy <- matrix(rnorm(10000), ncol=25)
#' set.seed(1000)
#' sing_vec <- getSVD(dummy, k = 1, sing.vec = "left")
#' 

getSVD <- function (matrix, k = 1, sing.vec = c("left", "right")) {  
  #center the matrix
  matrix <- t(scale(t(matrix), center = TRUE, scale = FALSE))
  #run SVD
  #set a random seed
  set.seed(1000)
  sing.vec <- match.arg(sing.vec)
  p.mat <- switch(sing.vec,
                  left=runSVD(matrix, k=k, BSPARAM=IrlbaParam())$u,
                  right=runSVD(matrix, k=k, BSPARAM=IrlbaParam())$v)
  #sum up the matrix
  csums <- colSums(matrix, na.rm=TRUE)
  #check for negative correlation
  #flip sign as necessary to ensure signal is associated
  #with pos. values
  if (cor(csums, p.mat) < 0) {
    p.mat <- -p.mat
  }
  #Chromosome length normalization
  p.mat <- p.mat * sqrt(length(p.mat))
  return(p.mat)
}
