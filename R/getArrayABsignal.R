#' @title Estimate A/B compartments from methylation array data
#'
#' @description 
#' \code{getArrayABsignal} returns estimated A/B compartments from methylation array data.
#'
#' @details This function is modified from the minfi::compartments to infer A/B compartments from array data
#' 
#'
#' @param obj Input GenomicRatioSet object 
#' @param res Compartment resolution (in bp)
#' @param parallel Should the inference be done in parallel?
#' @param allchrs Whether all autosomes should be used for A/B inference
#' @param chr Specific chromosomes to analyze
#' @param ... Additional arguments
#'
#' @return A p x n matrix (samples as columns and compartments as rows) of compartments
#' @import minfi
#' @export

getArrayABsignal <- function(obj, res=1e6, parallel=FALSE, allchrs=FALSE, chr = NULL, ...) {
  globalMeanSet <- .getMeanGrSet(obj)
  columns <- colnames(obj)
  names(columns) <- columns 
  
  getComp <- .getPairedArray
  if (allchrs == TRUE) getComp <- .getPairedAllChrsArray
  
  if (parallel) {
    options(mc.cores=detectCores()/2) # RAM blows up otherwise 
    do.call(cbind, 
            mclapply(columns,getComp,grSet=obj,globalMeanSet=globalMeanSet,chr=chr))
  } else { 
    do.call(cbind, 
            lapply(columns,getComp,grSet=obj,globalMeanSet=globalMeanSet,chr=chr))
  } 
}

#Helper function to get compartments as in minfi::compartments
.arraycompartments <- function(obj, resolution = 1e6, what="OpenSea",
                         chr = "chr22", method = c("pearson", "spearman"),
                         keep = FALSE) {
  stopifnot(length(chr) == 1 && chr %in% seqlevels(obj))
  method <- match.arg(method)
  gr <- .createCorMatrix(
    obj = obj,
    resolution = resolution,
    what = what,
    chr = chr,
    method = method)
  gr <- .extractAB(gr, keep = keep)
  gr$compartment <- .extractOpenClosedArray(gr)
  gr
}

#Helper function to perform mean imputation of missing values
.imputeMatrix <- function(matrix) {
    matrix[is.infinite(matrix) & matrix > 0] <- max(matrix[is.finite(matrix)])
    matrix[is.infinite(matrix) & matrix < 0] <- min(matrix[is.finite(matrix)])

    # Imputation of the missing values:
    # Probably should be using knn instead of mean imputation
    # Will leave for now but should be updated later
    missing <- which(is.na(matrix), arr.ind = TRUE)
    if (length(missing) != 0) {
        for (j in seq_len(nrow(missing))) {
            mean <- mean(
                x = matrix[missing[j, 1L], ][
                    is.finite(matrix[missing[j, 1L], ])],
                na.rm = TRUE)
            matrix[missing[j, 1L], missing[j, 2L]] <- mean
        }
    }
    matrix
}

#Helper function to return binned M-values
.returnBinnedMatrix <- function(gr.unbinnedCor, resolution){

    bin2D <- function(matrix, ids, n) {
        unique.ids <- sort(unique(ids))
        bin.matrix <- matrix(0, n, n)
        m <- length(unique.ids)
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                # TODO: We should be able to speed this one up a lot
                indices1 <- which(ids == unique.ids[i])
                indices2 <- which(ids == unique.ids[j])
                bin.matrix[unique.ids[i],unique.ids[j]] <- median(
                    x = matrix[indices1, indices2],
                    na.rm = TRUE)
            }
        }
        # TODO: This line should not be necessary
        bin.matrix[is.na(bin.matrix)] <- 1
        bin.matrix
    }

    # TODO: Use Bioconductor infrastructure to get chromosome lengths and
    #       don't directly call structure() to create object
    chr.lengths <- structure(
        c(249250621L, 243199373L, 198022430L, 191154276L, 180915260L,
          171115067L, 159138663L, 146364022L, 141213431L, 135534747L,
          135006516L, 133851895L, 115169878L, 107349540L, 102531392L,
          90354753L, 81195210L, 78077248L, 59128983L, 63025520L, 48129895L,
          51304566L, 155270560L, 59373566L),
        .Names = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                   "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                   "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22", "chrX", "chrY"))
    seqlengths(gr.unbinnedCor) <- chr.lengths[seqlevels(gr.unbinnedCor)]

    stopifnot(length(seqlevels(gr.unbinnedCor)) == 1 &&
                  !is.na(seqlengths(gr.unbinnedCor)))
    stopifnot("cor.matrix" %in% names(mcols(gr.unbinnedCor)))

    gr.binnedCor <- tileGenome(
        seqlengths = seqlengths(gr.unbinnedCor),
        tilewidth = resolution,
        cut.last.tile.in.chrom = TRUE)
    ids <- subjectHits(findOverlaps(gr.unbinnedCor, gr.binnedCor))
    gr.binnedCor$cor.matrix <- bin2D(
        matrix = gr.unbinnedCor$cor.matrix,
        ids = ids,
        n = length(gr.binnedCor))
    gr.binnedCor
}

#Helper function to remove empty bins
.removeBadBins <- function(gr) {
    good.bins <- which(!colAlls(gr$cor.matrix, value = 0))
    if (length(good.bins) < nrow(gr$cor.matrix)) {
        gr <- gr[good.bins]
        gr$cor.matrix <- gr$cor.matrix[, good.bins]
    }
    gr
}

#Helper function to extract whether the compartment is open or closed
.extractOpenClosedArray <- function(gr, cutoff = 0){
    pc <- gr$pc
    ifelse(pc < cutoff, "open", "closed")
}

#Helper function to get the first principal component
.getFirstPCarray <- function(matrix, method){
    # Centre the matrix
    center <- rowMeans2(matrix, na.rm = TRUE)
    matrix <- sweep(matrix, 1L, center, check.margin = FALSE)
    .fsvd(matrix, k = 1, method = method)$u
}

#Helper function to use a moving average smooth
.meanSmootherArray <- function(x, k = 1L, iter = 2L, na.rm = TRUE) {
    meanSmoother.internal <- function(x, k = 1L, na.rm = TRUE) {
        n <- length(x)
        y <- rep(NA_real_, n)

        window.mean <- function(x, j, k, na.rm = na.rm){
            if (k >= 1) {
                return(mean(x[seq(j - (k + 1L), j + k)], na.rm = na.rm))
            } else {
                x[j]
            }
        }

        for (i in (seq(k + 1L, n - k))) {
            y[i] <- window.mean(x, i, k, na.rm)
        }
        for (i in seq_len(k)) {
            y[i] <- window.mean(x, i, i - 1L, na.rm)
        }
        for (i in seq(n - k + 1L, n)) {
            y[i] <- window.mean(x, i, n - i, na.rm)
        }
        y
    }

    for (i in seq_len(iter)) {
        x <- meanSmoother.internal(x, k = k, na.rm = na.rm)
    }
    x
}

#Helper function to unitarize the A/B estimates
.unitarize <- function(x, medianCenter = TRUE) {
    if (medianCenter) x <- x - median(x, na.rm = TRUE)
    bad <- is.na(x)
    x[!bad] <- x[!bad] / sqrt(sum(x[!bad]^2))
    n.bad <- sum(bad)
    if (n.bad > 0) {
        message(
            sprintf("[.unitarize] %i missing values were ignored.\n", n.bad))
    }
    x
}

#Helper function to perform fsvd instead of the nipals method from mixOmics package
.fsvd <- function(A, k, i = 1, p = 2, method = c("qr", "svd", "exact")) {
    method <- match.arg(method)
    l <- k + p
    n <- ncol(A)
    m <- nrow(A)
    tall <- TRUE
    # NOTE: In the case a WIDE matrix is provided:
    if (m < n) {
        A <- t(A) # This is potentially expensive
        n <- ncol(A)
        m <- nrow(A)
        tall <- FALSE
    }
    if (l > ncol(A) && method != "exact") {
        stop("(k+p) is greater than the number of columns. Please decrease ",
             "the value of k.")
    }

    # Construct G to be n x l and Gaussian
    G <- matrix(rnorm(n * l, 0, 1), nrow = n, ncol = l)

    if (method == "svd") {
        # Power method:
        H <- A %*% G # m x l matrix
        for (j in seq_len(i)) {
            H <- A %*% (crossprod(A, H))
        }
        # NOTE: We use a SVD to find an othogonal basis Q:
        # TODO: Remove this commented line if not needed
        # H = FF %*% Omega %*% t(S)
        svd <- svd(crossprod(H))
        FF <- svd$u # l x l
        omega <- diag(1/sqrt(svd$d)) # l x l
        S <- H %*% FF %*% omega # m x l
        # Define the orthogonal basis:
        Q <- S[, seq_len(k), drop = FALSE] # m x k
        # TODO: Remove these commented lines if not needed
        # TT <- t(A) %*% Q # n x k
        # TT <- t(TT)
        TT <- crossprod(Q, A)
    } else if (method == "qr") {
        # NOTE: Need to create a list of H matrices
        h.list <- vector("list", i + 1L)
        h.list[[1]] <- A %*% G
        for (j in seq(2, i + 1L)) {
            h.list[[j]] <- A %*% (crossprod(A, h.list[[j - 1L]]))
        }
        H <- do.call("cbind", h.list) # n x [(1+1)l] matrix
        # QR algorithm
        Q <- qr.Q(qr(H, 0))
        # TODO: Remove these commented lines if not needed
        # TT <- t(A)%*%Q # n x [(i+1)l]
        # TT <- t(TT)
        TT <- crossprod(Q, A)
    }
    if (method == "svd" || method == "qr") {
        svd <- svd(TT)
        u <- Q %*% svd$u[,seq_len(k), drop = FALSE]
        v <- svd$v[, seq_len(k), drop = FALSE]
        d <- svd$d[seq_len(k)]
    } else {
        # Exact SVD
        svd <- svd(A)
        u <- svd$u[, seq_len(k), drop = FALSE]
        v <- svd$v[, seq_len(k), drop = FALSE]
        d <- svd$d[seq_len(k)]
    }
    if (!tall) {
        uu <- v
        v <- u
        u <- uu
    }
    list(u = u, v = v, d = d)
}

#Helper function to create a correlation matrix
.createCorMatrix <- function(obj, resolution = 100 * 1000, what = "OpenSea",
                            chr = "chr22", method = c("pearson", "spearman")) {
    stopifnot(length(chr) == 1 && chr %in% seqlevels(obj))
    method <- match.arg(method)

    if (is(obj, "GenomicMethylSet")) {
        obj <- ratioConvert(obj, what = "M", keepCN = FALSE)
    }
    assay(obj, "M") <- .imputeMatrix(getM(obj))

    # Next we subset to a chromosome, keep OpenSea probes and remove SNPs
    seqlevels(obj, pruning.mode = "coarse") <- chr
    obj <- obj[getIslandStatus(obj) %in% what,]
    obj <- dropLociWithSnps(obj, snps = c("CpG", "SBE"), maf = 0.01)

    gr.unbinnedCor <- granges(obj)
    gr.unbinnedCor$cor.matrix <- cor(t(getM(obj)), method = method)
    gr.cor <- .returnBinnedMatrix(gr.unbinnedCor, resolution = resolution)
    gr.cor <- .removeBadBins(gr.cor)
    gr.cor
}

#Helper function to extract A/B compartments
.extractAB <- function(gr, keep = TRUE, svdMethod = "qr"){
    if (!(is(gr, "GRanges") && "cor.matrix" %in% names(mcols(gr)))) {
        stop("'gr' must be an object created by createCorMatrix")
    }
    pc <- .getFirstPCarray(gr$cor.matrix, method = svdMethod)
    pc <- .meanSmootherArray(pc)
    pc <- .unitarize(pc)
    # Fix sign of eigenvector
    if (cor(colSums2(gr$cor.matrix), pc) < 0 ) {
        pc <- -pc
    }
    pc <- pc * sqrt(length(pc))
    gr$pc <- pc

    if (!keep) gr$cor.matrix <- NULL
    gr
}

.getMeanGrSet <- function(grSet) { 
  meanBeta <- matrix(rowMeans(getBeta(grSet), na.rm=TRUE), ncol=1) 
  subGrSet <- grSet[,1]
  assays(subGrSet)$Beta <- meanBeta
  sampleNames(subGrSet) <- "globalMean"
  return(subGrSet) 
}

.getPairedArray <- function(column, grSet, globalMeanSet=NULL, res=1e6, ...) {
  message("Computing shrunken compartment eigenscores for ", column, "...") 
  if(is.null(globalMeanSet)) globalMeanSet <- getMeanGrSet(grSet)
  .arraycompartments(cbind(grSet[,column], globalMeanSet), keep=FALSE, resolution=res, ...)$pc
}

.getPairedAllChrsArray <- function(column, grSet, globalMeanSet=NULL, res=1e6, ...) {
  if (is.null(globalMeanSet)) globalMeanSet <- getMeanGrSet(grSet)
  chrs <- intersect(paste0("chr", 1:22), seqlevels(grSet))
  names(chrs) <- chrs
  getPairedChr <- function(chr) { 
    message("Computing shrunken eigenscores for ", column, " on ", chr, "...") 
    .arraycompartments(cbind(grSet[,column],globalMeanSet),keep=FALSE,resolution=res,chr=chr)$pc
  }
  unlist(lapply(chrs, getPairedChr))
}
