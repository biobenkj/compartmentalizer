#' Filter cells using either k-means or Dirichlet process means clustering of sparsity metrics
#'
#' @name filterCells
#'
#' @param sparsity.mat A matrix of summarized sparsity measures
#' @param rse.obj The unfiltered RangedSummarizedExperiment object
#' @param cluster.method Clustering method to use (default: kmeans)
#' @param clusters How many clusters to generate; if NULL it will autopick the cluster number (default: NULL)
#' @param tol The tolerance or minimum difference in fraction of between cluster sum of squares over total for k-means auto-picking cluster number (default: 0.1)
#' @param plot.data Whether to plot the data
#' @param invert Invert which cluster is used to filter
#'
#' @return A filtered RangedSummarizedExperiment object and/or plot of the filtered data
#' @export
#' @import RColorBrewer
#' @import grid
#' @import SummarizedExperiment
#' @import ggplot2
#' @import viridis
#'
#' @examples
#'

filterCells <- function(sparsity.mat, rse.obj, cluster.method = c("kmeans", "dpmeans"),
                        clusters = NULL, tol = 0.1, plot.data = FALSE, invert = FALSE) {
  #this function is to filter out poor/low quality cells given the sparsity metrics
  #the underlying idea is that we will have 2 clusters of cells
  #1. cells with sufficient signal to use in our inference
  #2. cells with low/too sparse of signal genome-wide
  #input comes from getBinSparsity()
  #output will be a filtered RangeSummarizedExperiment object
  #ideally we can do this with lambda means - but not implemented yet...
  
  #check if the rse.obj is a RangedSummarizedExperiment object
  if (!is(rse.obj, "RangedSummarizedExperiment")) stop(rse.obj, " does not look like a RangedSummarizedExperiment object.")
  
  #get the clustering method
  cluster.method <- match.arg(cluster.method)
  
  #perform the clustering
  message("Clustering the data using: ", cluster.method)
  #the following only holds with kmeans
  #dp means does this for us and even better with lambda means - once that is implemented
  if (is.null(clusters)) {
    message("Auto-picking number of clusters based on sum of squares.")
    possible.clusters <- seq(2, 20)
    #initialize at 2 clusters
    cluster.choice <- 2
    #initialize ss for 1 cluster
    ss.choice <- 0
    for (c in possible.clusters) {
      #cluster with kmeans using the starting point
      kmeans.results <- stats::kmeans(t(sparsity.mat), centers = c)
      #compute sum of squares
      ss <- kmeans.results$betweenss / kmeans.results$totss
      message("Using ", c, " clusters: fraction between cluster sum of squares over total = ", ss*100, "%")
      if (ss.choice == 0) {
        ss.choice <- ss
        }
      else {
        ss.choice <- ifelse((ss - ss.choice) > tol, ss, "converged")
        #cluster.choice <- ifelse((ss - ss.choice) > 0.1, c, "converged")
      }
      if (as.character(ss.choice) == "converged") {
        message("Converged with clusters = ", c - 1)
        break()
      }
    }
    clusters <- c - 1
  }
  if (cluster.method == "kmeans") cluster.results <- stats::kmeans(t(sparsity.mat), centers = clusters)
  if (cluster.method == "dpmeans") stop("DPmeans is not implemented yet.")
  
  #filter the data
  if (!is.null(colnames(rse.obj))) {
    colname.match <- match(colnames(rse.obj), names(cluster.results$cluster))
    assay(rse.obj) <- assay(rse.obj)[,colname.match]
    }
  message("Filtering cells with too sparse of signal from ", deparse(substitute(rse.obj)), ".")

  #get the cluster to filter by sorting the per-sample sparsity measures and picking the lowest one
  #and querying which cluster that cell/sample belongs to
  cluster.to.filter <- cluster.results$cluster[names(sort(rowSums(t(sparsity.mat))))[1]]
  message("Filtering out cells belonging to cluster ", as.character(cluster.to.filter))
  if (invert) {
    filter.vec <- ifelse(cluster.results$cluster == cluster.to.filter, TRUE, FALSE) #make sure logic is inverted
  } else {
    filter.vec <- ifelse(cluster.results$cluster == cluster.to.filter, FALSE, TRUE) #make sure logic is inverted
  }
  filtered.data <- rse.obj[,filter.vec]
  message("Filtered out ", table(filter.vec)["FALSE"], " cells.")
  message("Kept ", table(filter.vec)["TRUE"], " cells with sufficient signal.")
  
  #plot the data
  if (plot.data) {
    message("Plotting results.")
    
    #make a boxplot 
    bxdata <- data.frame(cmeans = colMeans(sparsity.mat),
                         Clusters = factor(cluster.results$cluster))
    bplot <- ggplot(bxdata, aes(x = Clusters, y = cmeans, group = Clusters, color = Clusters)) +
      geom_boxplot() +
      geom_jitter() +
      scale_color_viridis(discrete = TRUE) +
      xlab("Clusters") +
      ylab("Mean genome-wide sparsity measures") +
      theme_bw(10) +
      ggtitle("Cluster assignment as a function of genome-wide sparsity measures") +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)
      )
    
    print(bplot)
    
    pca_sparsity <- prcomp(t(sparsity.mat))
    pr_comps <- data.frame(pca_sparsity$x)
    
    # Combine for plotting
    pr_comps$clusters <- factor(cluster.results$cluster)
    
    pca_plot <- ggplot(pr_comps, aes(x=PC1, y=PC2, color=clusters)) + 
      geom_point(size=3.5) + 
      ylim(-40, 40) +
      xlim(-40, 40) +
      scale_color_viridis(discrete = TRUE) + 
      theme_bw(10)
    
    # Plot percent variation explained
    prop_var <- data.frame(t(summary(pca_sparsity)$importance))
    names(prop_var) = c('sd', 'prop', 'cum')
    prop_var$num = 1:nrow(prop_var)
    
    var_plot <- ggplot(prop_var, aes(x=num, y=prop)) + 
      geom_point(size=3.5) + 
      geom_line() + 
      scale_x_continuous(limits = c(1, 22), breaks = 1:22) +
      xlab("Principal Component") + 
      ylab("Prop. of Variance") +
      theme_bw(10) +
      theme(
        axis.title.y = element_text(vjust=1),
        plot.margin = unit(c(0,0,0,6), "mm")
      )
    
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(4, 100)))
    print(pca_plot, vp = vplayout(1:3, 3:100))
    print(var_plot, vp = vplayout(4, 1:92))
  }
  
  return(filtered.data)
}
