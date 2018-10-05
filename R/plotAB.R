#' Plots A/B compartment estimates on a per chromosome basis
#'
#' Plot A/B compartments bins
#'
#' @param x      The GRanges object returned from getABSignal
#' @param main    Title for the plot
#' @param ylim      Y-axis limits (default is -1 to 1)
#' @param unitarize   Should the data be unitarized?
#' @param reverse    Reverse the sign of the PC values?
#' @param top.col    Top (pos. PC values) chromatin color to be plotted
#' @param bot.col    Bottom (neg. PC values) chromatin color to be plotted
#' 
#' @return    invisibly, the compartment estimates from the plot
#' 
#' @import    GenomicRanges
#' 
#' @export 
#' 
#' @examples 
#' 
#' library(GenomicRanges)
#' library(Homo.sapiens)
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
#' 
#' #Plot the A/B signal
#' par(mar=c(1,1,1,1))
#' par(mfrow=c(1,1))
#' plotAB(absignal$pc, ylim = c(-0.2, 0.2), unitarize = TRUE)
#' 
#' \dontrun{
#' #If plotting individual A/B signals using output from getCompartments
#' #Note: this function currently only supports plotting individual chromosomes from single samples
#' bin.chr1.ab <- getCompartments(data, "array", chrs = "chr1", genome = "hg19")
#' 
#' #For 7 samples
#' #Adjust ylim as necessary
#' par(mar=c(1,1,1,1))
#' par(mfrow=c(7,1))
#' plotAB(bin.chr1.ab[,1], ylim = c(-0.2, 0.2), unitarize = TRUE)
#' plotAB(bin.chr1.ab[,2], ylim = c(-0.2, 0.2), unitarize = TRUE, top.col = "goldenrod")
#' plotAB(bin.chr1.ab[,3], ylim = c(-0.2, 0.2), unitarize = TRUE, top.col = "darkblue")
#' plotAB(bin.chr1.ab[,4], ylim = c(-0.2, 0.2), unitarize = TRUE, top.col = "red")
#' plotAB(bin.chr1.ab[,5], ylim = c(-0.2, 0.2), unitarize = TRUE, top.col = "black")
#' plotAB(bin.chr1.ab[,6], ylim = c(-0.2, 0.2), unitarize = TRUE, top.col = "cyan")
#' plotAB(bin.chr1.ab[,7], ylim = c(-0.2, 0.2), unitarize = TRUE, top.col = "seagreen")
#' }



plotAB <- function(x, main="",ylim=c(-1, 1), unitarize=FALSE, reverse=FALSE,
                      top.col = "deeppink4", bot.col = "grey50"){
  if (unitarize){
    x <- .unitarize(x)
  }
  x <- as.numeric(x)
  if (reverse){
    x <- -x
  }
  
  n <- length(x)
  col <- rep(top.col, n)
  col[x<0] <- bot.col
  barplot(x, ylim=ylim, 
          bty="n", xlab="", ylab="",border=col, col=col, main=main)
}
