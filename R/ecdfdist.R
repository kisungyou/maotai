#' Distance Measures between Multiple Empirical Cumulative Distribution Functions
#' 
#' We measure distance between two empirical cumulative distribution functions (ECDF). For 
#' simplicity, we only take an input of \code{\link[stats]{ecdf}} objects from \pkg{stats} package. 
#' 
#' @param elist a length \eqn{N} list of \code{ecdf} objects.
#' @param method name of the distance/dissimilarity measure. Case insensitive (default: \code{ks}).
#' @param p exponent for \code{Lp} or \code{Wasserstein} distance. 
#' @param as.dist a logical; \code{TRUE} to return \code{dist} object, \code{FALSE} to return an \eqn{(N\times N)} symmetric matrix of pairwise distances (default: \code{FALSE}).
#' @param useR a logical; \code{TRUE} to use R implementation, \code{FALSE} to use C++ implementation (default: \code{TRUE}).
#' 
#' @return either \code{dist} object of an \eqn{(N\times N)} symmetric matrix of pairwise distances by \code{as.dist} argument.
#' 
#' @seealso \code{\link[stats]{ecdf}}
#' 
#' @examples
#' \donttest{
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
#' opar = par(no.readonly=TRUE)
#' par(pty="s")
#' image(dm[,nrow(dm):1], axes=FALSE, main=mks)
#' par(opar)
#' }
#' 
#' @export
ecdfdist <- function(elist, method=c("KS","Lp","Wasserstein"), p=2, as.dist=FALSE, useR=TRUE){
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
  if (as.logical(useR)){
    output = switch(mymethod, 
                    "ks"          = dist_ks(elist),
                    "wasserstein" = dist_wasserstein(elist, myp),
                    "lp"          = dist_lp(elist, myp))
  } else {
    output = switch(mymethod, 
                    "ks"          = fast_dist_ks(elist),
                    "wasserstein" = fast_dist_wasserstein(elist, myp),
                    "lp"          = fast_dist_lp(elist, myp))
  }
  
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

#' @keywords internal
#' @noRd
fast_dist_ks <- function(elist){
  trflist <- elist_fform(elist)        # returns list(tseq, fval=list of numeric)
  FF <- do.call(cbind, trflist$fval)   # L x N
  cpp_dist_ks(FF)
}

#' @keywords internal
#' @noRd
fast_dist_lp <- function(elist, p){
  trflist <- elist_fform(elist)
  tseq <- as.numeric(trflist$tseq)
  FF <- do.call(cbind, trflist$fval)
  cpp_dist_lp(tseq, FF, p)
}

#' @keywords internal
#' @noRd
fast_dist_wasserstein <- function(elist, p){
  nlist <- length(elist)
  qseq  <- base::seq(from = 1e-6, to = 1 - 1e-6, length.out = 8128L)
  # columns are quantiles for each ecdf on qseq
  Q <- vapply(elist,
              function(f) as.numeric(stats::quantile(f, qseq, names = FALSE)),
              numeric(length(qseq)))
  Q <- as.matrix(Q)                    # Lq x N
  cpp_dist_wasserstein(qseq, Q, p)
}









# ## toy example : random and uniform distributions
# n_test = 200
# n_half = round(n_test/2)
# mylist = vector("list", length=n_test)
# for (i in 1:n_half){
#   mylist[[i]] = stats::ecdf(stats::rnorm(500, sd=2))
# }
# for (i in (n_half+1):n_test){
#   mylist[[i]] = stats::ecdf(stats::runif(500, min=-5))
# }
# 
# # compute Kolmogorov-Smirnov distance
# microbenchmark::microbenchmark(
#   R_ks = ecdfdist(mylist, method="KS"),
#   R_lp = ecdfdist(mylist, method="Lp", p=2),
#   R_wt = ecdfdist(mylist, method="Wasserstein", p=2),
#   C_ks = ecdfdist(mylist, method="KS", useR=FALSE),
#   C_lp = ecdfdist(mylist, method="Lp", p=2, useR=FALSE),
#   C_wt = ecdfdist(mylist, method="Wasserstein", p=2, useR=FALSE),
#   times=3
# )
