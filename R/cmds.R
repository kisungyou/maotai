#' Classical Multidimensional Scaling
#' 
#' Classical multidimensional scaling aims at finding low-dimensional structure 
#' by preserving pairwise distances of data. 
#' 
#' @param data an \eqn{(n\times p)} matrix whose rows are observations.
#' @param ndim an integer-valued target dimension.
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
#' idata = as.matrix(iris[,1:4])
#' icol  = as.factor(iris[,5])   # class information
#' 
#' ## run Classical MDS
#' iris.cmds = cmds(idata, ndim=2)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' plot(iris.cmds$embed, col=icol, 
#'      main=paste0("STRESS=",round(iris.cmds$stress,4)))
#' par(opar)
#' 
#' @references 
#' \insertRef{torgerson_multidimensional_1952}{maotai} 
#' 
#' @export
cmds <- function(data, ndim=2){
  ############################################################
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* cmds : an input 'data' should be a matrix without any missing/infinite values.")
  }
  xdiss = stats::as.dist(cpp_pdist(data))
  
  ############################################################
  # Run and Return
  mydim  = round(ndim)
  output = hidden_cmds(xdiss, ndim=mydim)
  return(output)
}