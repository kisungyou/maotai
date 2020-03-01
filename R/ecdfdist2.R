#' Pairwise Measures for Two Sets of Empirical CDFs
#' 
#' We measure distance between two sets of empirical cumulative distribution functions (ECDF). For 
#' simplicity, we only take an input of \code{\link[stats]{ecdf}} objects from \pkg{stats} package. 
#' 
#' @param elist1 a length \eqn{M} list of \code{ecdf} objects.
#' @param elist2 a length \eqn{N} list of \code{ecdf} objects.
#' @param method name of the distance/dissimilarity measure. Case insensitive.
#' @param p exponent for \code{Lp} or \code{Wasserstein} distance. 
#' 
#' @return an \eqn{(M\times N)} matrix of pairwise distances.
#' 
#' @seealso \code{\link[stats]{ecdf}} \code{\link{ecdfdist}}
#' 
#' @examples 
#' ## toy example
#' #  first list : 10 of random and uniform distributions
#' mylist1 = list()
#' for (i in 1:10){ mylist1[[i]] = stats::ecdf(stats::rnorm(50, sd=2))}
#' for (i in 11:20){mylist1[[i]] = stats::ecdf(stats::runif(50, min=-5))}
#' 
#' #  second list : 15 uniform and random distributions
#' mylist2 = list()
#' for (i in 1:15){ mylist2[[i]] = stats::ecdf(stats::runif(50, min=-5))}
#' for (i in 16:30){mylist2[[i]] = stats::ecdf(stats::rnorm(50, sd=2))}
#' 
#' ## compute Kolmogorov-Smirnov distance
#' dm2ks = ecdfdist2(mylist1, mylist2, method="KS")
#' dm2lp = ecdfdist2(mylist1, mylist2, method="lp")
#' dm2wa = ecdfdist2(mylist1, mylist2, method="wasserstein")
#' nrs   = nrow(dm2ks)
#' 
#' ## visualize
#' opar = par(mfrow=c(1,3), no.readonly=TRUE)
#' image(dm2ks[,nrs:1], axes=FALSE, main="Kolmogorov-Smirnov")
#' image(dm2lp[,nrs:1], axes=FALSE, main="L2")
#' image(dm2wa[,nrs:1], axes=FALSE, main="Wasserstein")
#' par(opar)
#' 
#' @export
ecdfdist2 <- function(elist1, elist2, method=c("KS","Lp","Wasserstein"), p=2){
  ###############################################
  # Preprocessing
  if (!elist_check(elist1)){stop("* ecdfdist2 : input 'elist1' should be a list of 'ecdf' objects.")}
  if (!elist_check(elist2)){stop("* ecdfdist2 : input 'elist2' should be a list of 'ecdf' objects.")}
  methodss = c("ks","wasserstein","lp")
  mymethod = tolower(method)
  mymethod = match.arg(mymethod, methodss)
  myp      = as.integer(p)
  if (myp <= 0){
    stop("* ecdfdist2 : exponent 'p' should be a nonnegative number.")
  }
  
  ###############################################
  # Computation
  output = switch(mymethod, 
                  "ks"          = dist2_ks(elist1, elist2),
                  "wasserstein" = dist2_wasserstein(elist1, elist2, myp),
                  "lp"          = dist2_lp(elist1, elist2, myp))
  
  ###############################################
  # Return
  return(output)
}



# single functions --------------------------------------------------------
# (1) dist2_ks          : kolmogorov-smirnov
# (2) dist2_wasserstein : 1d wasserstein distance
# (3) dist2_lp          : Lp distance

#' @keywords internal
#' @noRd
dist2_ks <- function(elist1, elist2){
  M = length(elist1)
  N = length(elist2)
  
  trflst = elist_fform(c(elist1, elist2))
  flist1 = trflst$fval[1:M]
  flist2 = trflst$fval[(M+1):(M+N)]
  
  output = array(0,c(M,N))
  for (i in 1:M){
    fi = flist1[[i]]
    for (j in 1:N){
      fj = flist2[[j]]
      theval = max(abs(fi-fj))
      output[i,j] <- theval[1]
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
dist2_lp <- function(elist1, elist2, p){
  M = length(elist1)
  N = length(elist2)
  
  trflst = elist_fform(c(elist1, elist2))
  flist1 = trflst$fval[1:M]
  flist2 = trflst$fval[(M+1):(M+N)]
  
  output = array(0,c(M,N))
  for (i in 1:M){
    fi = flist1[[i]]
    for (j in 1:N){
      fj = flist2[[j]]
      if (is.infinite(p)){
        output[i,j] = base::max(base::abs(fi-fj))[1]
      } else {
        output[i,j] <- ((integrate_1d(trflst$tseq, (abs(fi-fj)^p)))^(1/p))
      }
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
dist2_wasserstein <- function(elist1, elist2, p){
  M = length(elist1)
  N = length(elist2)
  
  trflst = elist_fform(c(elist1, elist2))
  flist1 = trflst$fval[1:M]
  flist2 = trflst$fval[(M+1):(M+N)]
  
  qseq    = base::seq(from=1e-6, to=1-(1e-6), length.out=8128)
  quants1 = list()  # compute quantile functions first
  quants2 = list()
  for (i in 1:M){
    quants1[[i]] = as.double(stats::quantile(elist1[[i]], qseq))
  }
  for (j in 1:N){
    quants2[[j]] = as.double(stats::quantile(elist2[[j]], qseq))
  }
  
  output = array(0,c(M,N))
  for (i in 1:M){
    vali = quants1[[i]]
    for (j in 1:N){
      valj  = quants2[[j]]
      valij = abs(vali-valj)
      
      if (is.infinite(p)){
        output[i,j] = base::max(valij)
      } else {
        output[i,j] <- ((integrate_1d(qseq, valij^p))^(1/p))
      }
    }
  }
  return(output)
}