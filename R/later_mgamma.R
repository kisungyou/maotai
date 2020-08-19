#' Multivariate Gamma
#' 
#' \deqn{\Gamma_m (a)}
#' 
#' @keywords internal
#' @noRd
mgamma <- function(m, a, log=FALSE){
  m = round(m)
  if (length(a)==1){
    logval = (base::sum(base::lgamma(a - 0.5*((1:m)-1))))*(pi^(m*(m-1)/4))
    if (log){
      return(logval)
    } else {
      return(base::exp(logval))
    }
  } else {
    if (length(a)!=m){
      stop("* mgamma : for a vector-valued 'a', its length must be equal to 'm'.")
    }
    logval = base::exp(base::sum(base::lgamma(a - 0.5*((1:m)-1))))*(pi^(m*(m-1)/4))
    if (log){
      return(logval)
    } else {
      return(base::exp(logval))
    }
  }
}