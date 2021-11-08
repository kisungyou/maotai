#' Convert Covariance into Correlation Matrix
#' 
#' Given a covariance matrix, return a correlation matrix that has unit diagonals. 
#' We strictly impose (and check) whether the given input is a symmetric matrix 
#' of full-rank.
#' 
#' @param mat a \eqn{(p\times p)} covariance matrix.
#' 
#' @return a \eqn{(p\times p)} correlation matrix.
#' 
#' @examples 
#' \donttest{
#' # generate an empirical covariance scaled
#' prep_mat = stats::cov(matrix(rnorm(100*10),ncol=10))
#' prep_vec = diag(as.vector(stats::runif(10, max=5)))
#' prep_cov = prep_vec%*%prep_mat%*%prep_vec
#' 
#' # compute correlation matrix
#' prep_cor = cov2corr(prep_cov)
#' 
#' # visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(prep_cov, axes=FALSE, main="covariance")
#' image(prep_cor, axes=FALSE, main="correlation")
#' par(opar)
#' }
#' 
#' @export
cov2corr <- function(mat){
  # checker
  if (!check_covariance(mat)){
    stop("* cov2corr : an input 'mat' is not a valid covariance matrix.")
  }
  dvec = diag(1/sqrt(diag(mat)))
  return(dvec%*%mat%*%dvec)
}


# checker -----------------------------------------------------------------
#' @keywords internal
#' @noRd
check_covariance <- function(mat){
  # matrix
  if (!is.matrix(mat)){
    return(FALSE)
  }
  # symmetric
  if (!isSymmetric(mat)){
    return(FALSE)
  }
  # all positive diagonal
  if (any(diag(mat)<=0)){
    return(FALSE)
  }
  if (as.integer(Matrix::rankMatrix(mat)) < base::nrow(mat)){
    return(FALSE)
  }
  return(TRUE)
}
