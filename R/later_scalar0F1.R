#' Hypergeometric 0F1
#' 
#' 
#' @keywords internal
#' @noRd
scalar0F1 <- function(a, z, method=c("series")){
  # PREPARE
  mymethod = ifelse(missing(method),"series",
                    match.arg(tolower(method),
                              c("series")))
  # COMPUTE
  output = switch(mymethod,
                  series = scalar0F1.series(a, z))
  return(output)
}
#' @keywords internal
#' @noRd
scalar0F1.series <- function(a, z){
  no.stop = TRUE
  Mval    = 1
  n       = 0
  while (no.stop){
    n     = n+1
    M.now = exp(n*log(z) - sum(log((a + seq(from=0, to=(n-1), by=1)))) - base::lfactorial(n))
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
