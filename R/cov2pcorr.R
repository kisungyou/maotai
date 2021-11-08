#' Convert Covariance into Partial Correlation Matrix
#' 
#' Given a covariance matrix, return a partial correlation matrix that has unit diagonals. 
#' We strictly impose (and check) whether the given input is a symmetric matrix 
#' of full-rank.
#' 
#' @param mat a \eqn{(p\times p)} covariance matrix.
#' 
#' @return a \eqn{(p\times p)} partial correlation matrix.
#' 
#' @examples 
#' \donttest{
#' # generate an empirical covariance scaled
#' prep_mat = stats::cov(matrix(rnorm(100*10),ncol=10))
#' prep_vec = diag(as.vector(stats::runif(10, max=5)))
#' prep_cov = prep_vec%*%prep_mat%*%prep_vec
#' 
#' # compute correlation and partial correlation matrices
#' prep_cor = cov2corr(prep_cov)
#' prep_par = cov2pcorr(prep_cov)
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(prep_cov, axes=FALSE, main="covariance")
#' image(prep_cor, axes=FALSE, main="correlation")
#' image(prep_par, axes=FALSE, main="partial correlation")
#' par(opar)
#' }
#' 
#' @export
cov2pcorr <- function(mat){
  # checker
  if (!check_covariance(mat)){
    stop("* cov2corr : an input 'mat' is not a valid covariance matrix.")
  }
  return(src_cov2corr(mat))
}