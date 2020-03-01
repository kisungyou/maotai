#' Solve Lyapunov Equation
#' 
#' The Lyapunov equation is of form
#' \deqn{AX + XA^\top = Q}
#' where \eqn{A} and \eqn{Q} are square matrices of same size. Above form is also known as \emph{continuous} form. 
#' This is a wrapper of \code{armadillo}'s \code{sylvester} function.
#' 
#' @param A a \eqn{(p\times p)} matrix as above.
#' @param Q a \eqn{(p\times p)} matrix as above.
#' 
#' @return a solution matrix \eqn{X} of size \eqn{(p\times p)}.
#' 
#' @examples 
#' ## simulated example
#' #  generate square matrices
#' A = matrix(rnorm(25),nrow=5)
#' X = matrix(rnorm(25),nrow=5)
#' Q = A%*%X + X%*%t(A)
#' 
#' #  solve using 'lyapunov' function
#' solX = lyapunov(A,Q)
#' \dontrun{
#' pm1 = "* Experiment with Lyapunov Solver"
#' pm2 = paste("* Absolute Error  : ",norm(solX-X,"f"),sep="")
#' pm3 = paste("* Relative Error  : ",norm(solX-X,"f")/norm(X,"f"),sep="")
#' cat(paste(pm1,"\n",pm2,"\n",pm3,sep=""))
#' }
#' 
#' @references 
#' \insertRef{sanderson_armadillo_2016}{maotai}
#' 
#' \insertRef{eddelbuettel_rcpparmadillo_2014}{maotai}
#' 
#' @export
lyapunov <- function(A, Q){
  ###################################################################
  # check square matrix
  if (!check_sqmat(A)){
    stop("* lyapunov : an input 'A' should be a square matrix.")
  }
  if (!check_sqmat(Q)){
    stop("* lyapunov : an input 'Q' should be a square matrix.")
  }
  
  ###################################################################
  # arrange for RcppArmadillo format
  B = t(A)
  C = -Q
  
  ###################################################################
  # pass and return
  return(solve_lyapunov(A,B,C))
}