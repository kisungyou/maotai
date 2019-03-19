#' test for numerical approximation of gradient
#' 
#' 
#' @export
mygradF <- function(fn, x){
  return(gradF(fn,x))
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