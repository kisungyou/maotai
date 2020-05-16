#' Check for Distance Matrix
#' 
#' This function checks whether the distance matrix \eqn{D:=d_{ij} = d(x_i, x_j)} satisfies 
#' three axioms to make itself a semimetric, which are (1) \eqn{d_{ii} = 0}, (2) \eqn{d_{ij} > 0} for \eqn{i\neq j}, and 
#' (3) \eqn{d_{ij} = d_{ji}}.
#' 
#' @param d \code{"dist"} object or \eqn{(N\times N)} matrix of pairwise distances.
#' 
#' @return a logical; \code{TRUE} if it satisfies metric property, \code{FALSE} otherwise.
#' 
#' @examples 
#' ## Let's use L2 distance matrix of iris dataset
#' data(iris)
#' dx = as.matrix(stats::dist(iris[,1:4]))
#' 
#' # perturb d(i,j) 
#' dy = dx  
#' dy[1,2] <- dy[2,1] <- 10
#' 
#' # run the algorithm
#' checkdist(dx)
#' checkdist(dy)
#' 
#' @seealso \code{\link{checkmetric}}
#' @export
checkdist <- function(d){
  if (inherits(d, "dist")){
    d = as.matrix(d)
  } else {
    if (!is.matrix(d)){
      stop("* checkdist : input 'd' should be a matrix.")
    }
  }
  # 1. square matrix
  if (nrow(d)!=ncol(d)){
    message(" checkdist : input 'd' is not a square matrix.")
    return(FALSE)
  }
  # 2. zero diagonals
  if (any(diag(d)!=0)){
    message(" checkdist : input 'd' has non-zero diagonals.")
    return(FALSE)
  }
  # 3. all positive elements
  if (any(d < 0)){
    message(" checkdist : input 'd' has negative values.")
    return(FALSE)
  }
  # 4. symmetric
  if (!base::isSymmetric(d)){
    message(" checkdist : input 'd' is not symmetric.")
    return(FALSE)
  }
  
  return(TRUE)
}