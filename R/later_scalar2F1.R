#' Gauss Confluent Hypergeometric Function of Scalar Argument
#' 
#' 
#' @references 
#' \insertRef{butler_laplace_2002}{maotai}
#' 
#' @keywords internal
#' @noRd
scalar2F1 <- function(a, b, c, z, method=c("series","integral")){
  # PREPARE
  # if (abs(z) >= 1){
  #   stop("* scalar2F1 : '|z| < 1' is required.")
  # }
  mymethod = ifelse(missing(method),"series",
                    match.arg(tolower(method),
                              c("laplace","integral","series")))

  # COMPUTE
  output = switch(mymethod,
                  integral = scalar2F1.integral(a,b,c,z),
                  series   = scalar2F1.series(a,b,c,z),
                  laplace  = scalar2F1.laplace(a,b,c,z))
  return(output)
}
#' @keywords internal
#' @noRd
scalar2F1.integral <- function(a,b,c,z){
  # conditions not met
  
  # INTEGRATION
  func.int <- function(y){
    return((y^(a-1))*((1-y)^(c-a-1))*((1-(z*y))^(-b)))
  }
  myeps = 10*.Machine$double.eps
  term1 = stats::integrate(func.int, lower=(10*.Machine$double.eps), upper=1)$value
  term2 = 1/base::beta(a,c-a)
  return(term1*term2)
}
#' @keywords internal
#' @noRd
scalar2F1.series <- function(a,b,c,z){
  no.stop = TRUE
  Mval    = 1
  n       = 0
  while (no.stop){
    n     = n+1
    term1 = n*log(z) + sum(log((a + seq(from=0, to=(n-1), by=1)))) + sum(log((b + seq(from=0, to=(n-1), by=1))))
    term2 = sum(log((c + seq(from=0, to=(n-1), by=1)))) + base::lfactorial(n)
    Mnow  = exp(term1-term2)
    Mval  = Mval + Mnow
    if (abs(Mnow) < 1e-10){
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
scalar2F1.laplace <- function(a,b,c,x){
  tau  = (x*(b-a)) - c
  yhat = (2*a)/(sqrt((tau^2) - (4*a*x*(c-b))) - tau)
  
  r21  = ((yhat^2)/a) + (((1-yhat)^2)/(c-a)) - exp((log(b) + 2*log(x) + 2*log(yhat) + 2*log(1-yhat))-(2*log(1-(x*yhat)) + log(a) + log(c-a)))
  
  log1 = (c-0.5)*log(c)
  log2 = -0.5*log(r21)
  log3 = a*(log(yhat)-log(a)) + (c-a)*(log(1-yhat)-log(c-a)) + -b*log(1-(x*yhat))
  return(exp(log1+log2+log3))
}

# # special case
# myz = runif(1, min=-1, max=1)
# asin(myz)/myz
# scalar2F1(1/2,1/2,3/2,(myz^2), method = "series")
# scalar2F1(1/2,1/2,3/2,(myz^2), method = "integral")
# scalar2F1(1/2,1/2,3/2,(myz^2), method = "laplace")
