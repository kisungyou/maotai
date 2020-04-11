#' K-Means++ Clustering Algorithm
#' 
#' \eqn{k}-means++ algorithm is known to be a smart, careful initialization 
#' technique. It is originally intended to return a set of \eqn{k} points 
#' as initial centers though it can still be used as a rough clustering algorithm 
#' by assigning points to the nearest points.
#' 
#' @param data an \eqn{(n\times p)} matrix whose rows are observations.
#' @param k the number of clusters.
#' 
#' @return a length-\eqn{n} vector of class labels.
#' 
#' @examples 
#' ## use simple example of iris dataset
#' data(iris) 
#' mydata = as.matrix(iris[,1:4])
#' mycol  = as.factor(iris[,5])
#' 
#' ## find the low-dimensional embedding for visualization
#' my2d = cmds(mydata, ndim=2)$embed
#' 
#' ## apply 'kmeanspp' with different numbers of k's.
#' k2 = kmeanspp(mydata, k=2)
#' k3 = kmeanspp(mydata, k=3)
#' k4 = kmeanspp(mydata, k=4)
#' k5 = kmeanspp(mydata, k=5)
#' k6 = kmeanspp(mydata, k=6)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3))
#' plot(my2d, col=k2, main="k=2", pch=19, cex=0.5)
#' plot(my2d, col=k3, main="k=3", pch=19, cex=0.5)
#' plot(my2d, col=k4, main="k=4", pch=19, cex=0.5)
#' plot(my2d, col=k5, main="k=5", pch=19, cex=0.5)
#' plot(my2d, col=k6, main="k=6", pch=19, cex=0.5)
#' plot(my2d, col=mycol, main="true cluster", pch=19, cex=0.5)
#' par(opar)
#' 
#' @references 
#' \insertRef{arthur_kmeans_2007}{maotai}
#' 
#' @export
kmeanspp <- function(data, k=2){
  ############################################################
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* kmeanspp : an input 'data' should be a matrix without any missing/infinite values.")
  }
  xdiss = stats::as.dist(cpp_pdist(data))
  myk   = round(k)
  
  ############################################################
  # Run and Return
  output = hidden_kmeanspp(xdiss,k=myk)$cluster
  return(output)
}