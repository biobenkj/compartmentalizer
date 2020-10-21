#' Summarize the bootstrap compartment estimates and compute Agresti-Coull confidence intervals
#'
#' @name summarizeBootstraps
#'
#' @param boot.list List of bootstraps as GRanges objects 
#' @param est.ab The original compartment calls
#' @param q Confidence interval to compute (0.95 for 95 percent CI)
#' @param assay Type of assay we are working with
#'
#' @return A GRanges object with bootstraps summarized
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 

summarizeBootstraps <- function(boot.list, est.ab, q = 0.95,
                                assay = c("array", "atac", "bisulfite")) {
  #go through the estimated A/B compartments and compute proportions from the boot.list
  est.ab$score <- est.ab$pc
  message("Summarizing bootstraps.")
  #initialize open and closed counts to enumerate
  #est.ab$boot.open <- 0
  #est.ab$boot.closed <- 0
  #filter out failed compartment estimates
  boot.list <- removeEmptyBoots(boot.list)
  #summarize
  #create a dummy matrix and then rowSum them up
  boot.summary.mat.lst <- lapply(boot.list, function(b) {
    #add the pc to the granges object
    b$score <- b$pc
    #generate a dummy GRanges object
    est.ab.dummy <- est.ab
    est.ab.dummy$boot.open <- 0
    est.ab.dummy$boot.closed <- 0
    #convert to binary result for proportions
    #logic for methylation
    #open = eigen < 0
    #closed = eigen > 0
    #the logic is flipped for ATAC
    #open = eigen > 0
    #closed = eigen < 0
    if (assay == "atac") {
      b$open <- ifelse(b$score > 0, 1, 0)
      b$closed <- ifelse(b$score < 0, 1, 0)
    } else {
      b$open <- ifelse(b$score < 0, 1, 0)
      b$closed <- ifelse(b$score > 0, 1, 0)
    }
    #overlap by common intervals
    ol <- findOverlaps(b, est.ab.dummy)
    mcols(est.ab.dummy)$boot.open[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.open[subjectHits(ol)] + mcols(b)$open[queryHits(ol)]
    mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] <- mcols(est.ab.dummy)$boot.closed[subjectHits(ol)] + mcols(b)$closed[queryHits(ol)]
    #return the dummy mcols for bootstrapped open and closed calls
    return(as.matrix(cbind(mcols(est.ab.dummy)$boot.open,
                           mcols(est.ab.dummy)$boot.closed)))
  })
  
  #eunumerate bootstraps
  est.ab$boot.open <- rowSums(do.call("cbind", lapply(boot.summary.mat.lst, function(x) x[,1])))
  est.ab$boot.closed <- rowSums(do.call("cbind", lapply(boot.summary.mat.lst, function(x) x[,2])))
  
  message("Computing Agresti-Coull 95% confidence intervals.")
  est.ab$conf.est <- 0
  est.ab$conf.est.upperCI <- 0
  est.ab$conf.est.lowerCI <- 0
  conf.int <- lapply(1:length(est.ab), function(e) {
    #compute intervals
    if (est.ab[e,]$score > 0) {
      if (assay == "atac") {
        #assumes open for ATAC
        est <- agrestiCoullCI(est.ab[e,]$boot.open,
                              est.ab[e,]$boot.closed,
                              q = 0.95)
      } else {
        #assumes closed for arrays
        est <- agrestiCoullCI(est.ab[e,]$boot.closed,
                              est.ab[e,]$boot.open,
                              q = 0.95)
      }
      
      #this really isn't a good idea...
      #so... don't do it?
      #just return the ests
      #it will be the same length as est.ab
      return(est)
      #est.ab[e,]$conf.est.lowerCI <<- est[1]
      #est.ab[e,]$conf.est <<- est[2]
      #est.ab[e,]$conf.est.upperCI <<- est[3]
    }
    if (est.ab[e,]$score < 0) {
      if (assay == "atac") {
        #assumes closed for ATAC
        est <- agrestiCoullCI(est.ab[e,]$boot.closed,
                              est.ab[e,]$boot.open,
                              q = 0.95)
      } else {
        #assumes open for arrays
        est <- agrestiCoullCI(est.ab[e,]$boot.open,
                              est.ab[e,]$boot.closed,
                              q = 0.95)
      }
      #this really isn't a good idea...
      #est.ab[e,]$conf.est.lowerCI <<- est[1]
      #est.ab[e,]$conf.est <<- est[2]
      #est.ab[e,]$conf.est.upperCI <<- est[3]
      return(est)
    }
  })
  
  #combine the conf.est results into something sensible and bind with est.ab
  conf.int.ests <- do.call("rbind", conf.int)
  #should be of the form:
  #conf.est.lowerCI conf.est conf.est.upperCI
  est.ab$conf.est.lowerCI <- conf.int.ests[,1]
  est.ab$conf.est <- conf.int.ests[,2]
  est.ab$conf.est.upperCI <- conf.int.ests[,3]
  return(est.ab)
}
