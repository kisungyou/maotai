#' Metric Depth
#' 
#' Compute the metric depth proposed by \insertCite{geenens_2023_StatisticalDepthAbstract;textual}{maotai}, which is 
#' one generalization of statistical depth function onto the arbitrary metric space. Our implementation assumes that 
#' given the multivariate data it computes the (empirical) depth for all observations using under the Euclidean regime.
#' 
#' @param data an \eqn{(n\times p)} matrix whose rows are observations.
#' @return a length-\eqn{n} vector of empirical metric depth values.
#' 
#' @examples 
#' \dontrun{
#' ## use simple example of iris dataset
#' data(iris) 
#' X <- as.matrix(iris[,1:4])
#' y <- as.factor(iris[,5])
#' 
#' ## compute the metric depth
#' mdX <- metricdepth(X)
#' 
#' ## visualize
#' #  2-d embedding for plotting by MDS
#' X2d <- maotai::cmds(X, ndim=2)$embed
#' 
#' # get a color code for the metric depth
#' pal = colorRampPalette(c("yellow","red"))
#' 
#' # draw
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(X2d, pch=19, main="by class", xlab="", ylab="", col=y)
#' plot(X2d, pch=19, main="by depth", xlab="", ylab="", col=pal(150)[order(mdX)])
#' legend("bottomright", col=pal(2), pch=19, legend=round(range(mdX), 2))
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @export
metricdepth <- function(data){
  ## PREPROCESSING 
  if (!check_datamat(data)){
    stop("* metricdepth : an input 'data' should be a matrix without any missing/infinite values.")
  }
  xdiss = stats::as.dist(cpp_pdist(data))
  
  ## RUN AND RETURN
  output = hidden_metricdepth(xdiss)
  return(output)
}
