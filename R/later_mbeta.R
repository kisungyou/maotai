#' Multivariate Beta
#' 
#' 
#' @keywords internal
#' @noRd
mbeta <- function(m, a, b, log=FALSE){
  m      = round(m)
  logval = mgamma(m,a,log=TRUE) + mgamma(m,b,log=TRUE) - mgamma(m,(a+b),log=TRUE)
  if (log){
    return(logval)
  } else {
    return(base::exp(logval))
  }
}