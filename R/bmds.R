#' Bayesian Multidimensional Scaling
#' 
#' A Bayesian formulation of classical Multidimensional Scaling is presented.
#' Even though this method is based on MCMC sampling, we only return maximum a posterior (MAP) estimate
#' that maximizes the posterior distribution. Due to its nature without any special tuning,
#' increasing \code{mc.iter} requires much computation.
#' 
#' @param data an \eqn{(n\times p)} matrix whose rows are observations.
#' @param ndim an integer-valued target dimension.
#' @param par.a hyperparameter for conjugate prior on variance term, i.e., \eqn{\sigma^2 \sim IG(a,b)}. Note that \eqn{b} is chosen appropriately as in paper.
#' @param par.alpha hyperparameter for conjugate prior on diagonal term, i.e., \eqn{\lambda_j \sim IG(\alpha, \beta_j)}. Note that \eqn{\beta_j} is chosen appropriately as in paper.
#' @param par.step stepsize for random-walk, which is standard deviation of Gaussian proposal.
#' @param mc.iter the number of MCMC iterations.
#' @param verbose a logical; \code{TRUE} to show iterations, \code{FALSE} otherwise.
#' 
#' @return a named list containing
#' \describe{
#' \item{embed}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{stress}{discrepancy between embedded and origianl data as a measure of error.}
#' }
#' 
#' @examples
#' \donttest{
#' ## use simple example of iris dataset
#' data(iris) 
#' idata = as.matrix(iris[,1:4])
#' 
#' ## run Bayesian MDS
#' #  let's run 10 iterations only.
#' iris.cmds = cmds(idata, ndim=2)
#' iris.bmds = bmds(idata, ndim=2, mc.iter=5, par.step=(2.38^2)) 
#' 
#' ## extract coordinates and class information
#' cx = iris.cmds$embed # embedded coordinates of CMDS
#' bx = iris.bmds$embed #                         BMDS
#' icol = iris[,5]      # class information
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,1))
#' mc = paste0("CMDS with STRESS=",round(iris.cmds$stress,4))
#' mb = paste0("BMDS with STRESS=",round(iris.bmds$stress,4))
#' plot(cx, col=icol,pch=19,main=mc)
#' plot(bx, col=icol,pch=19,main=mb)
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{oh_bayesian_2001a}{maotai}
#' 
#' @export
bmds <- function(data, ndim=2, par.a=5, par.alpha=0.5, par.step=1, mc.iter=8128, verbose=TRUE){
  ############################################################
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* bmds : an input 'data' should be a matrix without any missing/infinite values.")
  }
  xdiss = stats::as.dist(cpp_pdist(data))
  
  ############################################################
  # Run the Hidden Function
  mydim   = round(ndim)
  mya     = as.double(par.a)
  myalpha = as.double(par.alpha)
  mystep  = as.double(par.step)
  myiter  = round(mc.iter)
  myshow  = as.logical(verbose)
  
  output = hidden_bmds(xdiss, ndim=mydim, 
                       par.a=mya, par.alpha=myalpha, 
                       par.step=mystep, mc.iter=myiter, verbose=myshow)
  
  ############################################################
  # Return the output
  return(output)
}