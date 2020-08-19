#' Confluent Hypergeometric Function of Matrix Argument
#' 
#' 
#' @references 
#' \insertRef{butler_laplace_2002}{maotai}
#' 
#' @keywords internal
#' @noRd
matrix1F1 <- function(a, b, Z, method=c("laplace")){
  # PREPARE
  if ((!isSymmetric(Z))||(!is.matrix(Z))){
    stop("* matrix1F1 : input 'Z' should be a symmetric matrix.")
  }
  mymethod = ifelse(missing(method),"laplace",
                    match.arg(tolower(method),
                              c("laplace")))
  # COMPUTATION
  output = switch(mymethod,
                  laplace = matrix1F1.laplace(a,b,Z))
  return(output)
}
#' @keywords internal
#' @noRd
matrix1F1.laplace <- function(a, b, X){ # checked in 1-dimension
  # Preliminary
  p = base::nrow(X)
  vec.xi = base::eigen(X)$values
  vec.yi = rep(0,p)
  for (i in 1:p){
    xi        = vec.xi[i]
    vec.yi[i] = (2*a)/(b-xi+sqrt(((xi-b)^2) + (4*a*xi)))
  }
  matR11 = array(0,c(p,p)) 
  for (i in 1:p){
    yi = vec.yi[i]
    for (j in 1:p){
      yj = vec.yi[j]
      matR11[i,j] = ((yi*yj)/a) + ((1-yi)*(1-yj)/(b-a))
    }
  }
  veclast = rep(0,p)
  for (i in 1:p){
    xi = vec.xi[i]
    yi = vec.yi[i]
    veclast[i] = base::exp(a*(log(yi)-log(a)) + (b-a)*(log(1-yi)-log(b-a)) + (xi*yi))
  }
  
  # Main
  log1 = ((b*p) - (p*(p+1)/4))*log(b)
  log2 = -0.5*base::sum(base::log(matR11))
  log3 = base::sum(base::log(veclast))
  return(base::exp(log1+log2+log3))
}