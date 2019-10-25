#' Distance Measures between Multiple Samples using Gaussian Mixture Models
#' 
#' Taking multiple observations (a sample) as a unit of analysis requires 
#' a measure of discrepancy between samples. \code{distgmm} fits finite 
#' Gaussian mixture models to each sample and use the fitted model as 
#' a representation of a single sample. A single model is selected via 
#' Bayesian Information Criterion (BIC). 
#' 
#' @param datalist a length \eqn{N} list of samples. All elements of the list should be of same type, either \code{vector} or \code{matrix} of same dimension (number of columns).
#' @param method name of the distance/dissimilarity measure.
#' @param maxk maximum number of clusters to be fitted using GMM.
#' @param as.dist a logical; \code{TRUE} to return \code{dist} object, \code{FALSE} to return an \eqn{(N\times N)} symmetric matrix of pairwise distances.
#' 
#' @return either \code{dist} object of an \eqn{(N\times N)} symmetric matrix of pairwise distances by \code{as.dist} argument.
#' 
#' @examples 
#' ## let's try two-dimensional data of 30 samples
#' ## single or mixture of two and three gaussian distributions.
#' dlist = list()
#' for (i in 1:10){
#'   dlist[[i]] = matrix(rnorm(120),ncol=2) 
#' }
#' for (i in 11:20){
#'   A = matrix(rnorm(60,mean=-4),ncol=2)
#'   B = matrix(rnorm(60,mean= 4),ncol=2)
#'   dlist[[i]] = rbind(A,B)
#' }
#' for (i in 21:30){
#'   A = matrix(rnorm(40,mean=-4),ncol=2)
#'   B = matrix(rnorm(40),ncol=2)
#'   C = matrix(rnorm(40,mean= 4),ncol=2)
#'   dlist[[i]] = rbind(A,B,C)
#' }
#' 
#' ## compute pairwise distances, expecting (3 x 3) block structure.
#' mm = distgmm(dlist, maxk=5)
#' 
#' ## visualize
#' opar <- par(pty="s")
#' image(mm[,nrow(mm):1], main="3-block pattern as expected")
#' par(opar)
#' 
#' @export
distgmm <- function(datalist, method=c("L2"), maxk=5, as.dist=FALSE){
  #######################################################
  # Preprocessing : checkers
  if (!check_datalist(datalist)){
    stop("* distgmm : an input should be a list containing samples of same dimension.")
  }
  maxk   = round(maxk)
  method = match.arg(method)
  nlist  = length(datalist)
  
  #######################################################
  # Compute : GMM
  list.gmm = list()
  if (is.vector(datalist[[1]])){
    vec.flag = TRUE
    for (n in 1:nlist){
     list.gmm[[n]] = mclust::Mclust(datalist[[n]], G=1:maxk, modelNames="V", verbose=FALSE)$parameters
    }  
  } else {
    vec.flag = FALSE
    for (n in 1:nlist){
      list.gmm[[n]] = mclust::Mclust(datalist[[n]], G=1:maxk, verbose=FALSE, modelNames="VVV")$parameters
    }
  }
  
  #######################################################
  # Compute : Pairwise Distance
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    objA = list.gmm[[i]]
    for (j in (i+1):nlist){
      objB = list.gmm[[j]]
      if (vec.flag){
        theval = switch(method, 
                        "L2" = distgmm_l2_1d(objA, objB))
        output[i,j] <- output[j,i] <- theval
      } else {
        theval = switch(method,
                        "L2" = distgmm_l2_nd(objA, objB))
        output[i,j] <- output[j,i] <- theval
      }
    }
  }
  
  #######################################################
  if (as.dist){
    return(stats::as.dist(output))
  } else {
    return(output)  
  }
}


# use Mclust 'parameters' object ------------------------------------------
#' @keywords internal
#' @noRd
distgmm_l2_1d <- function(objA, objB){
  weightA = as.vector(objA$pro)
  muA     = matrix(objA$mean, ncol=1)
  covA    = array(0,c(1,1,length(weightA)))
  for (i in 1:length(weightA)){
    covA[,,i] = objA$variance$sigmasq[i]
  }
  weightB = as.vector(objB$pro)
  muB     = matrix(objB$mean, ncol=1)
  covB    = array(0,c(1,1,length(weightB)))
  for (i in 1:length(weightB)){
    covB[,,i] = objB$variance$sigmasq[i]
  }
  
  ## run CPP (same for both 1d and nd cases)
  cpp.res = cpp_pairwise_L2(muA, muB, covA, covB)
  A = cpp.res$A
  B = cpp.res$B
  C = cpp.res$C
  
  ## matrix multiplication
  term1 = base::sum(as.vector(A%*%weightA)*weightA)
  term2 = base::sum(as.vector(B%*%weightB)*weightB)
  term3 = -2*base::sum(as.vector(C%*%weightB)*weightA)
  
  ## return distance/ L2 needs to be taken square root.
  return(base::sqrt(term1+term2+term3))
}
#' @keywords internal
#' @noRd
distgmm_l2_nd <- function(objA, objB){
  weightA = as.vector(objA$pro)
  muA     = t(objA$mean)
  covA    = objA$variance$sigma
  
  weightB = as.vector(objB$pro)
  muB     = t(objB$mean)
  covB    = objB$variance$sigma
  
  if (length(dim(covA)) < 3){
    tmpA = covA
    covA = array(0,c(ncol(muA),ncol(muA),1))
    covA[,,1] = as.matrix(tmpA)
  }
  if (length(dim(covB)) < 3){
    tmpB = covB
    covB = array(0,c(ncol(muB),ncol(muB),1))
    covB[,,1] = as.matrix(tmpB)
  }
  
  ## run CPP (same for both 1d and nd cases)
  cpp.res = cpp_pairwise_L2(muA, muB, covA, covB)
  A = cpp.res$A
  B = cpp.res$B
  C = cpp.res$C
  
  ## matrix multiplication
  term1 = base::sum(as.vector(A%*%weightA)*weightA)
  term2 = base::sum(as.vector(B%*%weightB)*weightB)
  term3 = -2*base::sum(as.vector(C%*%weightB)*weightA)
  
  ## return distance/ L2 needs to be taken square root.
  return(base::sqrt(term1+term2+term3))
}



# # personal experiment -----------------------------------------------------
# x = list()
# for (i in 1:20){
#   x[[i]] = matrix(rnorm(300*2),ncol=2)
# }
# for (i in 21:40){
#   x[[i]] = rbind(matrix(rnorm(150*2,mean=-4),ncol=2), matrix(rnorm(150*2,mean=4),ncol=2))
# }
# for (i in 41:60){
#   x[[i]] = rbind(matrix(rnorm(100*2,mean=-4),ncol=2), matrix(rnorm(100*2),ncol=2), matrix(rnorm(150*2,mean=4),ncol=2))
# }
# mm = distgmm(x, maxk=10)
# image(mm[,nrow(mm):1])

# bestgmm <- function(dat){
#   # belows are all automatically implemented in mclustBIC
#   # # run mclustBIC
#   # opt_gmm <- (mclust::mclustBIC(dat, G=1:9, verbose = FALSE))
#   # colgmm  <- colnames(opt_gmm)
#   # rowgmm  <- 1:9
#   # 
#   # # extract mclustBIC information
#   # mm      <- matrix(opt_gmm, nrow=nrow(opt_gmm))
#   # mm[is.na(mm)] = -Inf
#   # idmax   <- as.integer(which(mm == max(mm), arr.ind = TRUE)) # show the
#   # 
#   # nclust  <- rowgmm[idmax[1]]
#   # vartype <- colgmm[idmax[2]]
#   
#   # run Mclust
#   runobj <- mclust::Mclust(dat, G=1:10, verbose=FALSE)
#   
#   # run GMM with prespecified results
#   output  = list()
#   output$weight = runobj$parameters$pro
#   output$mu     = t(runobj$parameters$mean)
#   output$cov    = runobj$parameters$variance
#   return(output)
# }
