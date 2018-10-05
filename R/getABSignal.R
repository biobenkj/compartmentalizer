#' Calculate Pearson correlations of smoothed eigenvectors
#' 
#' This function is used to generate a list x to be passed to getABSignal
#'
#' @param x      A list object from getCorMatrix
#' @param k    Value of k for smoothing (default = 2)
#' @param iter    Number of iterations for moving average smoothing (default = 2)
#' @param squeeze    Whether squeezing was used (implies Fisher's Z transformation)
#' 
#' @return    A list x to pass to getABSignal
#' 
#' @import    GenomicRanges
#' @import    mixOmics
#' @import    Homo.sapiens
#' 
#' @export 
#' 
#' @examples 
#' 
#' library(GenomicRanges)
#' library(Homo.sapiens)
#' library(mixOmics)
#' 
#' #Generate random genomic intervals of 1-1000 bp on chr1-22
#' #Modified from https://www.biostars.org/p/225520/
#' random_genomic_int <- data.frame(chr = rep("chr14", 100))
#' random_genomic_int$start <- apply(random_genomic_int, 1, function(x) { round(runif(1, 0, seqlengths(Homo.sapiens)[x][[1]]), 0) })
#' random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
#' random_genomic_int$strand <- "*"
#' 
#' #Generate random counts
#' counts <- rnbinom(1000, 1.2, 0.4)
#' 
#' #Build random counts for 10 samples
#' count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
#' colnames(count.mat) <- paste0("sample_", seq(1:10))
#' 
#' #Bin counts
#' bin.counts <- getBinMatrix(count.mat, makeGRangesFromDataFrame(random_genomic_int), chr = "chr14", genome = "hg19")
#' 
#' #Calculate correlations
#' bin.cor.counts <- getCorMatrix(bin.counts)
#' 
#' #Get A/B signal
#' absignal <- getABSignal(bin.cor.counts)

getABSignal <- function(x, k = 2, iter = 2, squeeze = FALSE){
  message("Calculating eigenvectors...")
  pc <- .getFirstPC(x$binmat.cor)
  if (squeeze) pc <- ifisherZ(pc)
  message(paste0("Smoothing with a k of ", k, " for ", iter, " iterations..."))
  pc <- .meanSmoother(pc, k=k, iter=iter)
  message("Done smoothing...")
  gr <- x$gr
  gr$pc <- pc
  gr$compartments <- .extractOpenClosed(pc)
  return(gr)
}

#Internal function
#Code modified from https://github.com/Jfortin1/scATAC_Compartments
.getFirstPC <- function (matrix, ncomp = 1) {  
  matrix <- t(scale(t(matrix), center = TRUE, scale = FALSE))
  if (ncomp > 1) {
    p.mat <- nipals(matrix, ncomp = ncomp)$p
  }
  else {
    p.mat <- nipals(matrix, ncomp = ncomp)$p[, 1]
    csums <- colSums(matrix, na.rm=TRUE)
    if (cor(csums, p.mat) < 0){
      p.mat <- -p.mat
    }
  }
  p.mat <- p.mat * sqrt(length(p.mat)) #Chromosome length normalization
  return(p.mat)
}

#Internal function
#Code modified from minfi
# Author: Jean-Philippe Fortin
# May 6th 2015

.meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
  meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
    n <- length(x)
    y <- rep(NA,n)
    
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

#Internal function
#Code modified from minfi
# Author: Jean-Philippe Fortin
# May 6th 2015
.extractOpenClosed <- function(pc, cutoff = 0){
  ifelse(pc < cutoff, "open", "closed")
}
