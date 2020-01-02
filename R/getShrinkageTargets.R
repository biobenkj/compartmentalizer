#' Get the specified samples to shrink towards instead of the global mean
#'
#' @name getShrinkageTargets
#'
#' @param obj Input matrix
#' @param group Samples/colnames to use for targeted shrinkage
#'
#' @return A matrix composed of samples to shrink towards
#' @export
#'
#' @examples
#' dummy <- matrix(rnorm(1000), ncol=25)
#' dummy.sub <- getShrinkageTargets(dummy, group = c(1,5,8,10))
#' 

getShrinkageTargets <- function(obj, group) {
  if (is.null(colnames(obj))) {
    warning("No column names found to identify shrinkage targets.")
    warning("Assuming group values are column indices.")
    colnames(obj) <- seq_along(1:ncol(obj))
    }
  if (all(group %in% colnames(obj))) stargets.obj <- obj[,group]
  else (stop("Could not find ", group, " in the colnames of the input matrix..."))
  return(stargets.obj)
}