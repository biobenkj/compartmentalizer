#' Windowed mean smoother
#'
#' @name meanSmoother
#'
#' @param x Input data matrix where samples are columns and regions/loci are rows
#' @param k The number of windows to use (k=1 is 3 windows and k=2 is 5 windows)
#' @param iter The number of iterations to smooth
#' @param na.rm Whether to remove NAs prior to smoothing
#'
#' @return Smoothed data matrix
#' @export
#'
#' @examples
#' dummy <- matrix(rnorm(10000), ncol=25)
#' smooth.dummy <- meanSmoother(dummy)
#' 

meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
  meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
    if (k < 1) stop("k needs to be greater than or equal to 1...")
    if (length(x) < k) stop("Cannot smooth. Too few bins...")
    n <- length(x)
    y <- rep(NA,n)
    
    #note: k of 1 == 3 bins
    #note: k of 2 == 5 bins
    window.mean <- function(x, j, k, na.rm=na.rm){
      if (k>=1){
        return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
      } else {
        return(x[j])
      }    
    }
    
    for (i in (k+1):(n-k)){
      y[i] <- window.mean(x,i,k, na.rm)
    }
    for (i in 1:k){
      y[i] <- window.mean(x,i,i-1, na.rm)
    }
    for (i in (n-k+1):n){
      y[i] <- window.mean(x,i,n-i,na.rm)
    }
    y
  }
  
  for (i in 1:iter){
    x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
  }
  x
}
