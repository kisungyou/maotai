#' Solve Sylvester Equation
#' 
#' The Sylvester equation is of form
#' \deqn{AX + XB = C}
#' where \eqn{X} is the unknown and others are given. Though it's possible to have non-square \eqn{A} and \eqn{B} matrices, 
#' we currently support square matrices only. This is a wrapper of \code{armadillo}'s \code{sylvester} function.
#' 
#' @param A a \eqn{(p\times p)} matrix as above.
#' @param B a \eqn{(p\times p)} matrix as above.
#' @param C a \eqn{(p\times p)} matrix as above.
#' 
#' @return a solution matrix \eqn{X} of size \eqn{(p\times p)}.
#' 
#' @examples 
#' ## simulated example
#' #  generate square matrices
#' A = matrix(rnorm(25),nrow=5)
#' X = matrix(rnorm(25),nrow=5)
#' B = matrix(rnorm(25),nrow=5)
#' C = A%*%X + X%*%B
#' 
#' #  solve using 'sylvester' function
#' solX = sylvester(A,B,C)
#' pm1 = "* Experiment with Sylvester Solver"
#' pm2 = paste("* Absolute Error  : ",norm(solX-X,"f"),sep="")
#' pm3 = paste("* Relative Error  : ",norm(solX-X,"f")/norm(X,"f"),sep="")
#' cat(paste(pm1,"\n",pm2,"\n",pm3,sep=""))
#' 
#' 
#' @references 
#' \insertRef{sanderson_armadillo:_2016}{maotai}
#' 
#' \insertRef{eddelbuettel_rcpparmadillo:_2014}{maotai}
#' 
#' @export
sylvester <- function(A,B,C){
  ###################################################################
  # check square matrix
  if (!check_sqmat(A)){
    stop("* sylvester : an input 'A' should be a square matrix.")
  }
  if (!check_sqmat(B)){
    stop("* sylvester : an input 'B' should be a square matrix.")
  }
  if (!check_sqmat(C)){
    stop("* sylvester : an input 'C' should be a square matrix.")
  }
  
  ###################################################################
  # arrange and solve
  return(cpp_sylvester(A,B,-C))
}