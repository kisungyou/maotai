#' t-SNE Embedding
#' 
#' This function is a simple wrapper of \code{\link[Rtsne]{Rtsne}} function for 
#' t-Stochastic Neighbor Embedding for finding low-dimensional structure of 
#' the data embedded in the high-dimensional space. 
#' 
#' @param data an \eqn{(n\times p)} matrix whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param ... extra parameters to be used in \code{\link[Rtsne]{Rtsne}} function.
#' 
#' @return a named list containing
#' \describe{
#' \item{embed}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{stress}{discrepancy between embedded and origianl data as a measure of error.}
#' }
#' 
#' @examples 
#' ## use simple example of iris dataset 
#' data(iris) 
#' mydat = as.matrix(iris[,1:4])
#' mylab = as.factor(iris[,5])
#' 
#' ## run t-SNE and MDS for comparison
#' iris.cmds = cmds(mydat, ndim=2)
#' iris.tsne = tsne(mydat, ndim=2)
#' 
#' ## extract coordinates and class information
#' cx = iris.cmds$embed # embedded coordinates of CMDS
#' tx = iris.tsne$embed #                         t-SNE
#' 
#' ## visualize
#' #  main title
#' mc = paste("CMDS with STRESS=",round(iris.cmds$stress,4),sep="")
#' mt = paste("tSNE with STRESS=",round(iris.tsne$stress,4),sep="")
#' 
#' #  draw a figure
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,1))
#' plot(cx, col=mylab, pch=19, main=mc)
#' plot(tx, col=mylab, pch=19, main=mt)
#' par(opar)
#' 
#' @export
tsne <- function(data, ndim=2, ...){
  ############################################################
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* tsne : an input 'data' should be a matrix without any missing/infinite values.")
  }
  dx = (stats::as.dist(cpp_pdist(data)))
  kk = round(ndim)
  if ((length(kk)>1)||(kk<1)||(kk>=ncol(data))){
    stop("* tsne : 'ndim' should be an integer in [1,col(data)). ")
  }
  
  ############################################################
  # Run and Return
  output = hidden_tsne(dx, ndim=kk, ...)
  return(output)
}