#' Check for Metric Matrix
#' 
#' This function checks whether the distance matrix \eqn{D:=d_{ij} = d(x_i, x_j)} satisfies 
#' four axioms to make itself a semimetric, which are (1) \eqn{d_{ii} = 0}, (2) \eqn{d_{ij} > 0} for \eqn{i\neq j}, 
#' (3) \eqn{d_{ij} = d_{ji}}, and (4) \eqn{d_{ij} \leq d_{ik} + d_{kj}}.
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
#' checkmetric(dx)
#' checkmetric(dy)
#' 
#' @seealso \code{\link{checkdist}}
#' @export
checkmetric <- function(d){
  if (inherits(d, "dist")){
    d = as.matrix(d)
  } else {
    if (!is.matrix(d)){
      stop("* checkmetric : input 'd' should be a matrix.")
    }
  }
  # 1. square matrix
  if (nrow(d)!=ncol(d)){
    message(" checkmetric : input 'd' is not a square matrix.")
    return(FALSE)
  }
  # 2. zero diagonals
  if (any(diag(d)!=0)){
    message(" checkmetric : input 'd' has non-zero diagonals.")
    return(FALSE)
  }
  # 3. all positive elements
  if (any(d < 0)){
    message(" checkmetric : input 'd' contains negative values.")
    return(FALSE)
  }
  # 4. symmetric
  if (!base::isSymmetric(d)){
    message(" checkmetric : input 'd' is not symmetric.")
    return(FALSE)
  }
  # 5. triangle inequality
  return(cpp_triangle(d))
}

# data(iris)
# xx = as.matrix(iris[,1:4])
# dx = stats::dist(xx)
# dd = as.matrix(dx)
# 
# checkdist(dx)
# checkmetric(dx)
# 
# i=4
# j=11
# k=8
# 
# dd[i,j]
# dd[i,k]+dd[k,j]