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