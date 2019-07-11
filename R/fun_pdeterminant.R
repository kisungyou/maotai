#' Pseudo-Determinant
#' # (1) Holbrook's Note on Pseudo-Determinant -------------------------------
#     title : Differentiating The Pseudo Determinant
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
      del.new = del.old*multccc
      det.new = det(A+del.new*diag(n))/(del.new^(n-k))
      if ((abs(det.new-det.old)/abs(det.old))<1e-5){
        return(det.new)
      } else {
        del.old = del.new
        det.old = det.new
      }
    }
  }
}



# # personal tests ----------------------------------------------------------
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

