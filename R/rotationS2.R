#' Compute a Rotation on the 2-dimensional Sphere
#' 
#' A vector of unit norm is an element on the hypersphere. When two unit-norm 
#' vectors \eqn{x} and \eqn{y} in 3-dimensional space are given, this function 
#' computes a rotation matrix \eqn{Q} on the 2-dimensional sphere such that
#' \deqn{y=Qx}.
#' 
#' @param x a length-\eqn{3} vector. If \eqn{\|x\|\neq 1}, normalization is applied.
#' @param y a length-\eqn{3} vector. If \eqn{\|y\|\neq 1}, normalization is applied.
#' 
#' @return a \eqn{(3\times 3)} rotation matrix.
#' 
#' @examples
#' \donttest{
#' ## generate two data points
#' #  one randomly and another on the north pole
#' x = stats::rnorm(3)
#' x = x/sqrt(sum(x^2))
#' y = c(0,0,1)
#' 
#' ## compute the rotation
#' Q = rotationS2(x,y)
#' 
#' ## compare 
#' Qx = as.vector(Q%*%x)
#' 
#' ## print
#' printmat = rbind(Qx, y)
#' rownames(printmat) = c("rotated:", "target:")
#' print(printmat)
#' }
#' 
#' @export
rotationS2 <- function(x, y){
  ############################################################
  # Preprocessing 
  vec_x = as.vector(x); vec_x = vec_x/sqrt(sum(vec_x^2))
  vec_y = as.vector(y); vec_y = vec_y/sqrt(sum(vec_y^2))

  if (length(vec_x)!=3){    stop("rotationS2 : an input 'x' is not of length 3.")  }
  if (length(vec_y)!=3){    stop("rotationS2 : an input 'y' is not of length 3.")  }
  
  ############################################################
  # Run and Return
  output = rotateS2_main(vec_x, vec_y)
  return(output)
}



# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
rotateS2_main <- function(vec_u, vec_v){
  # orthonormal vector 'n'
  vec_n = pracma::cross(vec_u, vec_v)
  vec_n = vec_n/sqrt(sum(vec_n^2))
  
  # another orthonormal vector
  vec_t = pracma::cross(vec_n, vec_u)
  vec_t = vec_t/sqrt(sum(vec_t^2))
  
  # compute the angle
  alpha = base::atan2(sum(vec_v*vec_t), sum(vec_v*vec_u))
  
  # Rn(alpha)
  Rnalpha = cbind(c(base::cos(alpha), base::sin(alpha), 0),
                  c(-base::sin(alpha), base::cos(alpha), 0),
                  c(0,0,1))
  
  # T
  mat_T = cbind(vec_u,vec_t,vec_n)
  
  # rotator
  output = mat_T%*%Rnalpha%*%base::solve(mat_T)
  return(output)
}

