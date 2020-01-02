#' Filter to open sea CpG loci
#'
#' @name filterOpenSea
#'
#' @param obj Input SummarizedExperiment or GRanges object
#' @param genome Which genome to filter
#'
#' @return Filtered to open sea CpG loci
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' opensea <- filterOpenSea(array.data.chr14, genome = "hg19")
#' 

filterOpenSea <- function(obj, genome = c("hg19", "hg38", "mm10", "mm9"), other = NULL) {
  #get the desired open sea loci given the genome
  genome <- match.arg(genome)
  if (is.null(other)) {
    openseas.genome <- switch(genome,
                              hg19=data("openSeas.hg19", package="compartmap"),
                              hg38=data("openSeas.hg38", package="compartmap"),
                              mm10=data("openSeas.mm10", package="compartmap"),
                              mm9=data("openSeas.mm9", package="compartmap"))
  } else {
    #check if it's a GRanges flavored object
    if (!is(other, "GRanges")) stop("The 'other' input needs to be a GRanges of open sea regions")
    openseas.genome <- other
  }
  #Subset by overlaps
  message("Filtering to open sea CpG loci...")
  obj.openseas <- subsetByOverlaps(obj, get(openseas.genome))
  return(obj.openseas)
}

#' Gather open sea CpG from a GRanges of CpG islands
#' 
#' @description This function accepts a GRanges input of CpG islands that can
#' be derived from UCSC table browser and rtracklayer::import(yourbed.bed)
#'
#' @name filterOpenSea
#'
#' @param gr Input GRanges of CpG islands
#'
#' @return GRanges object that can be used with filterOpenSea()
#' @import rtracklayer
#' @import GenomicRanges
#' @export
#'
#' @examples
#' cpgi <- rtracklayer::import(system.file("inst/extdata/mm10_cpgi.bed", package = "compartmap"))
#' opensea_cpg <- getOpenSeas(cpgi)
#' 

getOpenSeas <- function(gr) {
  resorts <- trim(resize(gr, width(gr) + 8000, fix = "center"))
  openSeas <- subset(gaps(resorts), strand == "*")
  return(openSeas)
}

#' Retrieve probes, regions, etc. that fall within PMDs defined by Zhou et. al 2018
#' https://www.nature.com/articles/s41588-018-0073-4
#' 
#' @name getPMD
#'
#' @param obj Input SummarizedExperiment object
#' @param genome Which genome annotations to use
#'
#' @return A subsetted SummarizedExperiment containing loci within PMDs
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' pmd_loci <- getPMD(array.data.chr14, genome = "hg19")
#' 

getPMD <- function(obj, genome = c("hg19", "hg38")) {
  #make sure input is sane
  if (!checkAssayType(obj)) stop("Input needs to be a SummarizedExperiment.")
  #return PMDs derived from Zhou et al. 2018 Nature Genetics
  genome <- match.arg(genome)
  pmd.genome <- switch(genome,
                       hg19=data("PMDs.hg19", package = "compartmap"),
                       hg38=data("PMDs.hg38", package = "compartmap"))
  return(subsetByOverlaps(obj, pmd.genome))
}

#' Get the open and closed compartment calls based on sign of singular values
#'
#' @param gr Input GRanges with associated mcols that represent singular values
#' @param cutoff Threshold to define open and closed states
#' @param assay The type of assay we are working with
#'
#' @return A vector of binary/categorical compartment states
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' 
#' dummy <- matrix(rnorm(10000), ncol=25)
#' set.seed(1000)
#' sing_vec <- getSVD(dummy, k = 1, sing.vec = "left")
#' 
extractOpenClosed <- function(gr, cutoff = 0,
                              assay = c("array", "atac", "bisulfite")){
  #check for input to be GRanges
  if (!is(gr, "GRanges")) stop("Input needs to be a GRanges.")
  if (!("pc" %in% names(mcols(gr)))) stop("Need to have an mcols column be named 'pc'.")
  assay <- match.arg(assay)
  if (assay %in% c("array", "bisulfite")) return(ifelse(gr$pc < cutoff, "open", "closed"))
  if (assay == "atac") return(ifelse(gr$pc < cutoff, "closed", "open"))
}

#' Check if the assay is a SummarizedExperiment
#'
#' @param obj Input object
#'
#' @return Boolean
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' checkAssayType(array.data.chr14)

checkAssayType <- function(obj) {
  #helper function to check the class of an object
  is(obj, "SummarizedExperiment")
}

#' Get the assay names from a SummarizedExperiment object
#'
#' @param se Input SummarizedExperiment object
#'
#' @return The names of the assays
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' getAssayNames(array.data.chr14)

getAssayNames <- function(se) {
  #helper function to check the assay slot names
  names(assays(se))
}

#' Remove rows with NAs exceeding a threshold
#'
#' @param se Input SummarizedExperiment object
#' @param rowmax The maximum NAs allowed in a row as a fraction
#' @param assay The type of assay we are working with
#'
#' @return A filtered matrix
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' cleanAssayRows(array.data.chr14, assay = "array")

cleanAssayRows <- function(se, rowmax = 0.5,
                           assay = c("array", "atac", "bisulfite")) {
  assay <- match.arg(assay)
  switch(assay,
         array = se[rowMeans(is.na(assays(se)$Beta)) < rowmax,],
         atac = se[rowMeans(is.na(assays(se)$counts)) < rowmax,],
         bisulfite = se[rowMeans(is.na(assays(se)$counts)) < rowmax,])
}

#' Remove columns/cells/samples with NAs exceeding a threshold
#'
#' @param se Input SummarizedExperiment object
#' @param colmax The maximum number of NAs allowed as a fraction
#' @param assay The type of assay we are working with
#'
#' @return A filtered matrix
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' cleanAssayCols(array.data.chr.14, assay = "array")
cleanAssayCols <- function(se, colmax = 0.8,
                           assay = c("array", "atac", "bisulfite")) {
  assay <- match.arg(assay)
  switch(assay,
         array = se[,colMeans(is.na(assays(se)$Beta)) < colmax],
         atac = se[,colMeans(is.na(assays(se)$counts)) < colmax],
         bisulfite = se[,colMeans(is.na(assays(se)$counts)) < colmax])
}

#' Helper function: squeezed logit
#'
#' @param x       a vector of values between 0 and 1 inclusive
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between -Inf and +Inf
#'
#' @import        gtools
#'
#' @export 
flogit <- function(x, sqz=0.000001) {
  x[ which(x < sqz) ] <- sqz 
  x[ which(x > (1 - sqz)) ] <- (1 - sqz)
  logit(x)
}

#' Helper function: expanded expit
#'
#' @param x       a vector of values between -Inf and +Inf
#' @param sqz     the amount by which to 'squeeze', default is .000001
#'
#' @return        a vector of values between 0 and 1 inclusive
#'
#' @import        gtools
#'
#' @export 
fexpit <- function(x, sqz=0.000001) {
  (((((inv.logit(x) * 2) - 1) / (1 - sqz)) + 1) / 2)
}

#' Get the chromosomes from an object
#'
#' @param obj Input SummarizedExperiment object
#'
#' @return A character vector of chromosomes present in an object
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' data("meth_array_450k_chr14", package = "compartmap")
#' getChrs(array.data.chr14)
#' 

getChrs <- function(obj) {
  #get the chromosomes present in the object
  return(unique(as.character(seqnames(obj))))
}

#' Remove bootstrap estimates that failed
#'
#' @param obj Input list object with elements 'pc' and 'gr'
#'
#' @return A filtered list object
#' @export

removeEmptyBoots <- function(obj) {
  #remove NAs from a bootstrap list
  #this can happen if the correlation between the bins and eigenvector fails
  #theoretically we can recover these but need an additional utility to find consensus
  na.filt <- unlist(lapply(obj, function(n) ifelse(any(is.na(n)), FALSE, TRUE)))
  obj <- obj[na.filt]
  return(obj)
}
