#' Invert, or "fix", compartments that have a minimum confidence score (1-min.conf)
#'
#' @name fixCompartments
#'
#' @param obj Input RaggedExperiment or output of condenseSE
#' @param min.conf Minimum confidence score to use
#' @param parallel Whether to run in parallel
#' @param cores How many cores to use if running in parallel
#'
#' @return A "fixed" set of compartments
#' @import RaggedExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#'  

fixCompartments <- function(obj, min.conf = 0.8,
                            parallel = FALSE, cores = 1) {
  #this function will invert or "fix" compartments based on bootstrapping results
  if (is(obj, "RaggedExperiment")) obj <- condenseSE(obj, sample.name = colnames(assay(obj)))
  
  #if we somehow only have 1 sample
  if (is(obj, "GRanges")) {
    if (!any((names(mcols(obj)) %in% "conf.est"))) stop("Bootstrapping was not performed. Cannot fix compartments.")
    message("Assuming we only have a single sample to process.")
    invert_compartments <- apply(mcols(obj), 1, function(c) {
      return(ifelse(c["conf.est"] < 1-min.conf, TRUE, FALSE))
    })
    message("Fixing compartments using a minimum confidence score of ", min.conf*100, "%")
    mcols(obj)$flip.compartment <- invert_compartments
    #add a new column for flipped scores
    mcols(obj)$flip.score <- mcols(obj)$score
    #flip the score
    mcols(obj)$flip.score[invert_compartments] <- -(mcols(obj)$score[invert_compartments])
    #add a new column for flipped CIs
    mcols(obj)$flip.conf.est <- mcols(obj)$conf.est
    mcols(obj)$flip.conf.est.upperCI <- mcols(obj)$conf.est.upperCI
    mcols(obj)$flip.conf.est.lowerCI <- mcols(obj)$conf.est.lowerCI
    #flip the conf.est
    mcols(obj)$flip.conf.est[invert_compartments] <- 1-(mcols(obj)$conf.est[invert_compartments])
    #flip the upperCI
    mcols(obj)$flip.conf.est.upperCI[invert_compartments] <- 1-(mcols(obj)$conf.est.lowerCI[invert_compartments])
    mcols(obj)$flip.conf.est.lowerCI[invert_compartments] <- 1-(mcols(obj)$conf.est.upperCI[invert_compartments])
    return(obj)
  }
  
  message("Fixing compartments using a minimum confidence score of ", min.conf*100, "%")
  #go through and invert compartments based on the min.conf
  if (parallel) {
    flip_compartments_lst <- mclapply(obj, function(s) {
      if (!any((names(mcols(s)) %in% "conf.est"))) stop("Bootstrapping was not performed. Cannot fix compartments.")
      invert_compartments <- apply(mcols(s), 1, function(c) {
        return(ifelse(c["conf.est"] < 1-min.conf, TRUE, FALSE))
      })
      mcols(s)$flip.compartment <- invert_compartments
      #add a new column for flipped scores
      mcols(s)$flip.score <- mcols(s)$score
      #flip the score
      mcols(s)$flip.score[invert_compartments] <- -(mcols(s)$score[invert_compartments])
      #add a new column for flipped CIs
      mcols(s)$flip.conf.est <- mcols(s)$conf.est
      mcols(s)$flip.conf.est.upperCI <- mcols(s)$conf.est.upperCI
      mcols(s)$flip.conf.est.lowerCI <- mcols(s)$conf.est.lowerCI
      #flip the conf.est
      mcols(s)$flip.conf.est[invert_compartments] <- 1-(mcols(s)$conf.est[invert_compartments])
      #flip the upperCI
      mcols(s)$flip.conf.est.upperCI[invert_compartments] <- 1-(mcols(s)$conf.est.lowerCI[invert_compartments])
      mcols(s)$flip.conf.est.lowerCI[invert_compartments] <- 1-(mcols(s)$conf.est.upperCI[invert_compartments])
      return(s)
    }, mc.cores = cores)
  } else {
    flip_compartments_lst <- lapply(obj, function(s) {
      if (!any((names(mcols(s)) %in% "conf.est"))) stop("Bootstrapping was not performed. Cannot fix compartments.")
      invert_compartments <- apply(mcols(s), 1, function(c) {
        return(ifelse(c["conf.est"] < 1-min.conf, TRUE, FALSE))
      })
      mcols(s)$flip.compartment <- invert_compartments
      #add a new column for flipped scores
      mcols(s)$flip.score <- mcols(s)$score
      #flip the score
      mcols(s)$flip.score[invert_compartments] <- -(mcols(s)$score[invert_compartments])
      #add a new column for flipped CIs
      mcols(s)$flip.conf.est <- mcols(s)$conf.est
      mcols(s)$flip.conf.est.upperCI <- mcols(s)$conf.est.upperCI
      mcols(s)$flip.conf.est.lowerCI <- mcols(s)$conf.est.lowerCI
      #flip the conf.est
      mcols(s)$flip.conf.est[invert_compartments] <- 1-(mcols(s)$conf.est[invert_compartments])
      #flip the upperCI
      mcols(s)$flip.conf.est.upperCI[invert_compartments] <- 1-(mcols(s)$conf.est.lowerCI[invert_compartments])
      mcols(s)$flip.conf.est.lowerCI[invert_compartments] <- 1-(mcols(s)$conf.est.upperCI[invert_compartments])
      return(s)
    })
  }
  names(flip_compartments_lst) <- names(obj)
  return(flip_compartments_lst)
}
