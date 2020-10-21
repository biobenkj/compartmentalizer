#' Condense a RaggedExperiment to a list of SummarizedExperiments
#'
#' @param obj Input RaggedExperiment
#'
#' @return A list of SummarizedExperiments corresponding to the assays in the input
#' @import RaggedExperiment
#' @export
#'
#' @examples
#' grl <- GRangesList(GRanges(c("A:1-5", "A:4-6", "A:10-15"), score=1:3),
#' GRanges(c("A:1-5", "B:1-3"), score=4:5))
#' names(grl) <- c("A", "B")                  
#' x <- RaggedExperiment(grl)
#' x.condense <- condenseRE(x)
condenseRE <- function(obj) {
  #this is a function to extract relevant information from a RaggedExperiment
  #it will build a list of SummarizedExperiments with relevant information
  #from computing compartments
  if (!is(obj, "RaggedExperiment")) stop("Input needs to be a RaggedExperiment")
  se_list <- lapply(1:length(assayNames(obj)), function(a) {
    compactSummarizedExperiment(obj, i = a)
  })
  #do NOT use getAssayNames here
  #for some reason it causes memory to skyrocket
  names(se_list) <- assayNames(obj)
  return(se_list)
}

#' Condense the output of condenseRE to reconstruct per-sample GRanges objects to plot
#'
#' @param obj Output of condenseRE or can be a RaggedExperiment
#' @param sample.name Vector of samples/cells to extract
#'
#' @return GRanges or list of per-sample GRanges to pass to plotAB or export
#' @import RaggedExperiment
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' grl <- GRangesList(GRanges(c("A:1-5", "A:4-6", "A:10-15"), score=1:3),
#' GRanges(c("A:1-5", "B:1-3"), score=4:5))
#' names(grl) <- c("A", "B")                  
#' x <- RaggedExperiment(grl)
#' condense.x <- condenseSE(x, sample.name = "A")
condenseSE <- function(obj, sample.name = NULL) {
  if (is.null(sample.name)) stop("Need a sample/cell to extract via sample.name.")
  #condense the input to something that can be plotted with plotAB
  if (is(obj, "RaggedExperiment")) obj <- condenseRE(obj)
  #make sure there are some assays to work with
  if (length(obj) < 1) stop("No assays found to condense.")
  #check and make sure that the names needed are found in the column names
  if (!(all(sample.name %in% colnames(assay(obj[[1]]))))) stop("The sample.name(s) not found in the colnames of the assays.")
  #check and see how many samples we are extracting
  if (length(sample.name) == 1) {
    obj.dense <- lapply(1:length(obj), function(a) {
      gr.sub <- granges(obj[[a]])
      mcols(gr.sub) <- assay(obj[[a]])[,sample.name]
      names(mcols(gr.sub)) <- assayNames(obj[[a]])
      return(gr.sub)
    })
    return(Reduce("merge", obj.dense))
  } else {
    obj.dense.lst <- lapply(sample.name, function(s) {
      obj.dense <- lapply(1:length(obj), function(a) {
        gr.sub <- granges(obj[[a]])
        mcols(gr.sub) <- assay(obj[[a]])[,s]
        names(mcols(gr.sub)) <- assayNames(obj[[a]])
        return(gr.sub)
      })
      return(Reduce("merge", obj.dense))
    })
    names(obj.dense.lst) <- sample.name
    return(obj.dense.lst)
  }
}