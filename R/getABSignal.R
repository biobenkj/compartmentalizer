#' Calculate Pearson correlations of smoothed eigenvectors
#' 
#' This function is used to generate a list x to be passed to getABSignal
#'
#' @param x      A list object from getCorMatrix
#' @param squeeze    Whether squeezing was used (implies Fisher's Z transformation)
#' @param assay What kind of assay are we working on ("array", "atac", "bisulfite")
#' 
#' @return    A list x to pass to getABSignal
#' 
#' @import    SummarizedExperiment
#' 
#' @export 
#' 
#' @examples 
#' 
#' library(SummarizedExperiment)
#' library(Homo.sapiens)
#' library(BiocSingular)
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

getABSignal <- function(x, squeeze = FALSE,
                        assay = c("array", "atac", "bisulfite")){
  message("Calculating eigenvectors.")
  pc <- getSVD(x$binmat.cor, sing.vec = "right")
  if (squeeze) pc <- ifisherZ(pc)
  message("Smoothing eigenvector.")
  pc <- switch(assay,
               array = meanSmoother(pc, k=1, iter=2),
               atac = meanSmoother(pc, k=2, iter=2),
               bisulfite = meanSmoother(pc, k=1, iter=2))
  message("Done smoothing.")
  gr <- x$gr
  gr$pc <- pc
  gr$compartments <- extractOpenClosed(gr, assay = assay)
  return(gr)
}
