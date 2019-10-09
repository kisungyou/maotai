#' Numerical Approximation to Gradient of a Function with Matrix Argument
#' 
#' For a given function \eqn{f:\mathbf{R}^{n\times p} \rightarrow \mathbf{R}}, 
#' we use finite difference scheme that approximates a gradient at a given point \eqn{x}. 
#' In Riemannian optimization, this can be used as a proxy for 
#' ambient gradient. Use with care since it may accumulate numerical error.
#' 
#' @param fn a function that takes a matrix of size \eqn{(n\times p)} and returns a scalar value.
#' @param x an \eqn{(n\times p)} matrix where the gradient is to be computed.
#' @param h step size for centered difference scheme. 
#' 
#' @return an approximate numerical gradient matrix of size \eqn{(n\times p)}.
#' 
#' @examples 
#' ## function f(X) = <a,Xb> for two vectors 'a' and 'b'
#' #  derivative w.r.t X is ab'
#' #  take an example of (5x5) symmetric positive definite matrix
#' 
#' #  problem settings
#' a   <- rnorm(5)
#' b   <- rnorm(5)
#' ftn <- function(X){
#'   return(sum(as.vector(X%*%b)*a))
#' }       # function to be taken derivative
#' myX <- matrix(rnorm(25),nrow=5)  # point where derivative is evaluated
#' myX <- myX%*%t(myX)
#' 
#' # main computation
#' sol.true <- base::outer(a,b)
#' sol.num1 <- matderiv(ftn, myX, h=1e-1) # step size : 1e-1
#' sol.num2 <- matderiv(ftn, myX, h=1e-5) #             1e-3
#' sol.num3 <- matderiv(ftn, myX, h=1e-9) #             1e-5
#' 
#' ## visualize/print the results
#' expar = par(mfrow=c(2,2),pty="s")
#' image(sol.true, main="true solution")
#' image(sol.num1, main="h=1e-1")
#' image(sol.num2, main="h=1e-5")
#' image(sol.num3, main="h=1e-9")
#' par(expar)
#' 
#' ntrue = norm(sol.true,"f")
#' cat('* Relative Errors in Frobenius Norm ')
#' cat(paste("*  h=1e-1   : ",norm(sol.true-sol.num1,"f")/ntrue,sep=""))
#' cat(paste("*  h=1e-5   : ",norm(sol.true-sol.num2,"f")/ntrue,sep=""))
#' cat(paste("*  h=1e-9   : ",norm(sol.true-sol.num3,"f")/ntrue,sep=""))
#' 
#' @references 
#' \insertRef{kincaid_numerical_2009}{maotai}
#' 
#' @export
matderiv <- function(fn, x, h=0.001){
  if (h <= 0){
    stop("* matderiv : 'h' should be a nonnegative real number.")
  }
  hval = max(sqrt(.Machine$double.eps), h)
  return(gradF(fn,x,hval))
}
# h = 0.001
# X = matrix(rnorm(9),nrow=3)
# X = X%*%t(X)
# dX = array(0,c(3,3))
# fX = function(x){return(sum(diag(x%*%x)))}
# for (i in 1:3){
#   for (j in 1:3){
#     Xp = X
#     Xm = X
#     Xp[i,j] = Xp[i,j] + h
#     Xm[i,j] = Xm[i,j] - h
#     dX[i,j] = (fX(Xp)-fX(Xm))/(2*h)
#   }
# }
# dX