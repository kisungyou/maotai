#' Calculate the Pseudo-Determinant of a Matrix
#' 
#' When a given square matrix \eqn{A} is rank deficient, determinant is zero. 
#' Still, we can compute the pseudo-determinant by multiplying all non-zero 
#' eigenvalues. Since thresholding to determine near-zero eigenvalues is subjective, 
#' we implemented the function as of original limit problem. When matrix is 
#' non-singular, it coincides with traditional determinant.
#' 
#' @param A a square matrix whose pseudo-determinant be computed.
#' 
#' @return a scalar value for computed pseudo-determinant.
#' 
#' @examples 
#' ## show the convergence of pseudo-determinant
#' #  settings
#' n = 10
#' A = cov(matrix(rnorm(5*n),ncol=n))   # (n x n) matrix
#' k = as.double(Matrix::rankMatrix(A)) # rank of A
#' 
#' # iterative computation
#' ntry = 11
#' del.vec = exp(-(1:ntry))
#' det.vec = rep(0,ntry)
#' for (i in 1:ntry){
#'   del = del.vec[i]
#'   det.vec[i] = det(A+del*diag(n))/(del^(n-k))
#' }
#' 
#' # visualize the results
#' plot(1:ntry, det.vec, main=paste("true rank is ",k," out of ",n,sep=""),"b", xlab="iterations")
#' abline(h=pdeterminant(A),col="red",lwd=1.2)
#' 
#' 
#' @references 
#' \insertRef{holbrook_differentiating_2018}{maotai}
#' 
#' @export
pdeterminant <- function(A){ ## wrapped-up function
  if (!is.matrix(A)){
    stop("* pdeterminant : input 'A' should be a matrix.")
  }
  if (nrow(A)!=ncol(A)){
    stop("* pdeterminant : input 'A' should be a square matrix.")
  }
  n = nrow(A)
  k = as.double(Matrix::rankMatrix(A))
  if (k==n){
    return(base::determinant(A))
  } else {
    multccc = 0.9
    del.old = 1
    det.old = det(A+del.old*diag(n))/(del.old^(n-k))
    for (i in 1:496){
      # print(paste("iteration ",i," initiated...",sep=""))
      del.new = del.old*multccc
      det.new = det(A+del.new*diag(n))/(del.new^(n-k))
      if ((abs(det.new-det.old)/abs(det.old))<1e-5){
        return(det.new)
      } else {
        del.old = del.new
        det.old = det.new
      }
    }
    return(det.old)
  }
}



# # # personal tests ----------------------------------------------------------
# n = 10
# A = cov(matrix(rnorm(5*n),ncol=n))
# k = as.double(Matrix::rankMatrix(A))
# 
# ntry = 20
# del.vec = exp(-(1:ntry))
# det.vec = rep(0,ntry)
# for (i in 1:ntry){
#   del = del.vec[i]
#   det.vec[i] = det(A+del*diag(n))/(del^(n-k))
# }
# plot(1:ntry, det.vec, main=paste("true rank is ",k,"/",n,sep=""),"b")
# abline(h=pdeterminant(A),col="red",lwd=1.2)

