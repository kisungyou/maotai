#' Distance Measures between Multiple Empirical Cumulative Distribution Functions
#' 
#' We measure distance between two empirical cumulative distribution functions (ECDF). For 
#' simplicity, we only take an input of \code{\link[stats]{ecdf}} objects from \pkg{stats} package. 
#' 
#' @param elist a length \eqn{N} list of \code{ecdf} objects.
#' @param method name of the distance/dissimilarity measure. Case insensitive.
#' @param p exponent for \code{Lp} or \code{Wasserstein} distance. 
#' @param as.dist a logical; \code{TRUE} to return \code{dist} object, \code{FALSE} to return an \eqn{(N\times N)} symmetric matrix of pairwise distances.
#' 
#' @return either \code{dist} object of an \eqn{(N\times N)} symmetric matrix of pairwise distances by \code{as.dist} argument.
#' 
#' @seealso \code{\link[stats]{ecdf}}
#' 
#' @examples 
#' ## toy example : 10 of random and uniform distributions
#' mylist = list()
#' for (i in 1:10){
#'   mylist[[i]] = stats::ecdf(stats::rnorm(50, sd=2))
#' }
#' for (i in 11:20){
#'   mylist[[i]] = stats::ecdf(stats::runif(50, min=-5))
#' }
#' 
#' ## compute Kolmogorov-Smirnov distance
#' dm = ecdfdist(mylist, method="KS")
#' 
#' ## visualize
#' mks  =" KS distances of 2 Types"
#' opar = par(pty="s", no.readonly=TRUE)
#' image(dm[,nrow(dm):1], axes=FALSE, main=mks)
#' par(opar)
#' 
#' @export
ecdfdist <- function(elist, method=c("KS","Lp","Wasserstein"), p=2, as.dist=FALSE){
  ###############################################
  # Preprocessing
  if (!elist_check(elist)){
    stop("* ecdfdist : input 'elist' should be a list of 'ecdf' objects.")
  }
  methodss = c("ks","wasserstein","lp")
  mymethod = tolower(method)
  mymethod = match.arg(mymethod, methodss)
  myp      = round(p)
  if (myp <= 0){
    stop("* ecdfdist : exponent 'p' should be a nonnegative number.")
  }
  
  ###############################################
  # Computation
  output = switch(mymethod, 
                  "ks"          = dist_ks(elist),
                  "wasserstein" = dist_wasserstein(elist, myp),
                  "lp"          = dist_lp(elist, myp))
  
  ###############################################
  # Report
  if (as.dist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}



# single functions --------------------------------------------------------
# (1) dist_ks          : kolmogorov-smirnov
# (2) dist_wasserstein : 1d wasserstein distance
# (3) dist_lp          : Lp distance


#' @keywords internal
#' @noRd
dist_ks <- function(elist){
  trflist  = elist_fform(elist)
  flist = trflist$fval
  nlist = length(flist)
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    fi = flist[[i]]
    for (j in (i+1):nlist){
      fj = flist[[j]]
      theval = max(abs(fi-fj))
      output[i,j] <- output[j,i] <- theval[1]
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
dist_lp <- function(elist, p){
  nlist = length(elist)
  trflist  = elist_fform(elist)
  flist = trflist$fval
  nlist = length(flist)
  output = array(0,c(nlist,nlist))
  if (is.infinite(p)){
    for (i in 1:(nlist-1)){
      fi = flist[[i]]
      for (j in (i+1):nlist){
        fj = flist[[j]]
        output[i,j] <- output[j,i] <- base::max(base::abs(fi-fj))[1]
      }
    } 
  } else {
    for (i in 1:(nlist-1)){
      fi = flist[[i]]
      for (j in (i+1):nlist){
        fj = flist[[j]]
        theval = ((integrate_1d(trflist$tseq, (abs(fi-fj)^p)))^(1/p))
        output[i,j] <- output[j,i] <- theval
      }
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
dist_wasserstein <- function(elist, p){
  nlist = length(elist)
  qseq  = base::seq(from=1e-6, to=1-(1e-6), length.out=8128)
  quants = list() # compute quantile functions first
  for (i in 1:nlist){
    quants[[i]] = as.double(stats::quantile(elist[[i]], qseq))
  }
  
  
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    vali = quants[[i]]
    for (j in (i+1):nlist){
      valj = quants[[j]]
      valij = abs(vali-valj)
      
      if (is.infinite(p)){
        output[i,j] <- output[j,i] <- base::max(valij)
      } else {
        theval <- ((integrate_1d(qseq, valij^p))^(1/p))
        output[i,j] <- output[j,i] <- theval
      }
    }
  } 
  
  return(output)
}

## wasserstein : http://www-users.math.umn.edu/~bobko001/preprints/2016_BL_Order.statistics_Revised.version.pdf