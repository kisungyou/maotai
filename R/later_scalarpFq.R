#' General Form of Hypergeometric Function
#' 
#' 
#' @keywords internal
#' @noRd
scalarpFq <- function(veca, vecb, z){
  p = length(veca)
  q = length(vecb)
  
  no.stop = TRUE
  Mval    = 1
  n       = 0
  while (no.stop){
    n  = n+1
    terma = 0
    for (i in 1:p){
      terma = terma + sum(log((veca[i] + seq(from=0, to=(n-1), by=1))))
    }
    termb = 0
    for (j in 1:q){
      termb = termb + sum(log((vecb[j] + seq(from=0,to=(n-1),by=1))))
    }
    Mnow = exp(n*log(z) + terma - termb - base::lfactorial(n))
    Mval = Mnow + Mval
    if (abs(Mnow) < 1e-10){
      no.stop=FALSE
    } 
    if (n>100){
      no.stop=FALSE
    }
  }
  return(Mval)
}