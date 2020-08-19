#' Gauss Confluent Hypergeometric Function of Matrix Argument
#' 
#' @references 
#' \insertRef{butler_laplace_2002}{maotai}
#' 
#' @keywords internal
#' @noRd
matrix2F1 <- function(a, b, c, Z, method=c("laplace")){
  # PREPARE
  if ((!isSymmetric(Z))||(!is.matrix(Z))){
    stop("* matrix2F1 : input 'Z' should be a symmetric matrix.")
  }
  mymethod = ifelse(missing(method),"laplace",
                    match.arg(tolower(method),
                              c("laplace")))
  
  # RUN
  output = switch(mymethod,
                  laplace = matrix2F1.laplace(a,b,c,Z))
  return(output)
}
#' @keywords internal
#' @noRd
matrix2F1.laplace <- function(a, b, c, X){
  p      = base::nrow(X)
  vec.xx = base::eigen(X)$values
  vec.yy = rep(0,p)
  for (i in 1:p){
    tau       = (vec.xx[i]*(b-a)) - c
    vec.yy[i] = (2*a)/(sqrt((tau^2) - (4*a*vec.xx[i]*(c-b))) - tau)
  }
  vec.Li  = rep(0,p)
  for (i in 1:p){
    vec.Li[i] = (vec.xx[i]*vec.yy[i]*(1-vec.yy[i]))/(1-(vec.xx[i]*vec.yy[i]))
  }
  matR21 = array(0,c(p,p))
  for (i in 1:p){
    xi = vec.xx[i]
    yi = vec.yy[i]
    for (j in 1:p){
      xj = vec.xx[j]
      yj = vec.yy[j]
      
      term1 = exp(log(yi)+log(yj)-log(a))
      term2 = exp(log(1-yi)+log(1-yj)-log(c-a))
      term3top = log(b)+log(xi)+log(xj)+log(yi)+log(yj)+log(1-yi)+log(1-yj)
      term3bot = log(1-(xi*yi))+log(1-(xj*yj))+log(a)+log(c-a)
      term3 = exp(term3top-term3bot)
      matR21[i,j] = term1+term2-term3
    }
  }
  
  
  
  log1 = ((c*p) - (p*(p+1)/4))*log(c)
  log2 = -0.5*sum(log(matR21))
  log3 = (a*(log(vec.yy)-log(a)))+((c-a)*(log(1-vec.yy)-log(c-a)))-(b*log(1-(vec.xx*vec.yy)))
  return(base::exp(log1+log2+log3))
}



# # # special case
# myz = runif(1, min=-1, max=1)
# 
# asin(myz)/myz
# scalar2F1(1/2,1/2,3/2,(myz^2), method = "series")
# scalar2F1(1/2,1/2,3/2,(myz^2), method = "integral")
# scalar2F1(1/2,1/2,3/2,(myz^2), method = "laplace")
# 
# matrix2F1(1/2,1/2,3/2,matrix((myz^2)), method = "laplace")