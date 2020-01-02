# turn a matrix with appropriate rownames back into a GRanges
#
# note: the granges() method is already documented and this brings no change
#
setAs("matrix", "GRanges",
      function(from) { 
        gr <- as(rownames(from), "GRanges")
        mcols(gr) <- as(from, "DataFrame")
        return(gr) 
      })

# wrap
setMethod("granges", "matrix", function(x) as(x, "GRanges"))

# turn it back 
setAs("GRanges", "matrix", 
      function(from) { 
        mat <- as(mcols(from), "matrix")
        rownames(mat) <- as.character(from) 
        return(mat)  
      }) 

# wrap
setMethod("as.matrix", "GRanges", function(x, ...) as(x, "matrix"))
