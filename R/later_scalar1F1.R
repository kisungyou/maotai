#' Kummer's Confluent Hypergeometric Function of Scalar Argument
#' 
#' @references 
#' \insertRef{butler_laplace_2002}{maotai}
#' 
#' @keywords internal
#' @noRd
scalar1F1 <- function(a, b, z, method=c("series","laplace","integral")){
  # PREPARE
  mymethod = ifelse(missing(method),"series",
                    match.arg(tolower(method),
                              c("laplace","integral","series")))
  
  # COMPUTE
  output = switch(mymethod,
                  integral = scalar1F1.integral(a,b,z),
                  series   = scalar1F1.series(a,b,z),
                  laplace  = scalar1F1.laplace(a,b,z))
  return(output)
}
#' @keywords internal
#' @noRd
scalar1F1.integral <- function(a, b, z){
  # REQUIREMENT
  if (!((b>a)&&(a>0))){
    stop("scalar1F1 : 'integral' method requires 'b > a > 0'.")
  }
  # INTEGRATION
  func.int <- function(y){
    return((y^(a-1))*((1-y)^(b-a-1))*exp(z*y))
  }
  myeps = 10*.Machine$double.eps
  term1 = stats::integrate(func.int, lower=(10*.Machine$double.eps), upper=1)$value
  term2 = 1/base::beta(a,b-a)
  return(term1*term2)
}
#' @keywords internal
#' @noRd
scalar1F1.series <- function(a, b, z){
  no.stop = TRUE
  Mval    = 1
  n       = 0
  while (no.stop){
    n     = n+1
    M.now = exp(sum(log((a + seq(from=0, to=(n-1), by=1)))) + n*log(z) - sum(log((b + seq(from=0, to=(n-1), by=1)))) - base::lfactorial(n))
    Mval  = Mval + M.now
    if (abs(M.now) < 1e-10){
      no.stop = FALSE
    }
    if (n > 100){
      no.stop = FALSE
    }
  }
  return(Mval)
}
#' @keywords internal
#' @noRd
scalar1F1.laplace <- function(a, b, x){
  yhat = (2*a)/(b-x+sqrt(((x-b)^2) + (4*a*x)))
  r11  = (yhat^2)/a + ((1-yhat)^2)/(b-a)
  
  log1 = (b-0.5)*log(b)
  log2 = -0.5*log(r11)
  log3 = a*(log(yhat)-log(a))
  log4 = (b-a)*(log(1-yhat)-log(b-a))
  log5 = x*yhat
  
  output = exp(log1+log2+log3+log4+log5)
  return(output)
}

# mya = 1/2
# myb = sample(2:20, 1)/2
# myx = abs(1+rnorm(1, sd=10))
# 
# fAsianOptions::kummerM(myx,mya,myb)
# scalar1F1(mya,myb,myx,method="integral")
# scalar1F1(mya,myb,myx,method="series")
# scalar1F1(mya,myb,myx,method="laplace")
