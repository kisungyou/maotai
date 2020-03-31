#' Weiszfeld Algorithm for Computing L1-median
#' 
#' Geometric median, also known as L1-median, is a solution to the following problem
#' \deqn{\textrm{argmin} \sum_{i=1}^n \| x_i - y \|_2 }
#' for a given data \eqn{x_1,x_2,\ldots,x_n \in R^p}. 
#' 
#' @param X an \eqn{(n\times p)} matrix for \eqn{p}-dimensional signal. If vector is given, it is assumed that \eqn{p=1}.
#' @param weights \code{NULL} for equal weight \code{rep(1/n,n)}; otherwise, it has to be a vector of length \eqn{n}.
#' @param maxiter maximum number of iterations.
#' @param abstol stopping criterion 
#' 
#' @examples 
#' ## generate sin(x) data with noise for 100 replicates
#' set.seed(496)
#' t = seq(from=0,to=10,length.out=20)
#' X = array(0,c(100,20))
#' for (i in 1:100){
#'    X[i,] = sin(t) + stats::rnorm(20, sd=0.5)
#' }
#' 
#' ## compute L1-median and L2-mean
#' vecL2 = base::colMeans(X)
#' vecL1 = weiszfeld(X)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' matplot(t(X[1:5,]), type="l", main="5 generated data", ylim=c(-2,2))
#' plot(t, vecL2, type="l", col="blue", main="L2-mean",   ylim=c(-2,2))
#' plot(t, vecL1, type="l", col="red",  main="L1-median", ylim=c(-2,2))
#' par(opar)
#' 
#' @export
weiszfeld <- function(X, weights=NULL, maxiter=496, abstol=1e-6){
  ###############################################
  # Preprocessing
  if (is.vector(X)){
    X = matrix(X, ncol = 1)
  }
  n = nrow(X)
  d = ncol(X)
  
  if ((length(weights)==0)&&(is.null(weights))){
    weights = rep(1/n, n)
  }
  if ((!is.vector(weights))||(length(weights)!=n)){
    stop(paste0("* weiszfeld : 'weights' should be a vector of length ",n))
  }
  
  ###############################################
  # Prepare
  myiter = round(maxiter)
  mytol  = as.double(abstol)
  xinit  = as.vector(base::colMeans(X))
  epsnum = (100*.Machine$double.eps)
  
  ###############################################
  # Compute and Return
  return(as.vector(cpp_weiszfeld(X, mytol, myiter, xinit, weights, epsnum)))
}



# ## example from Gmedian:weiszfeld
# ## Robustness of the geometric median of n=3 points in dimension d=2.
# library(Gmedian)
# 
# ## Computation speed
# ## Simulated data - Brownian paths
# n <- 1e4
# d <- 50
# x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
# x <- t(apply(x,1,cumsum))
# 
# library(microbenchmark)
# microbenchmark::microbenchmark(
#   "resG1" = as.vector(Weiszfeld(x)$median),
#   "resG2" = as.vector(Gmedian(x)),
#   "resWW" = as.vector(maotai::weiszfeld(x))  
# )

