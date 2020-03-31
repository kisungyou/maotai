#' Generate Index for Moving Block Bootstrapping
#' 
#' Assuming data being dependent with cardinality \code{N}, \code{boot.mblock} returns 
#' a vector of index that is used for moving block bootstrapping.
#' 
#' @param N the number of observations.
#' @param b the size of a block to be drawn.
#' 
#' @return a vector of length \code{N} for moving block bootstrap sampling.
#' 
#' @examples 
#' ## example : bootstrap confidence interval of mean and variances
#' vec.x = seq(from=0,to=10,length.out=1000)
#' vec.y = sin(1.21*vec.x) + 2*cos(3.14*vec.x) + rnorm(1000,sd=1.5)
#' data.mu  = mean(vec.y)
#' data.var = var(vec.y)
#' 
#' ## apply moving block bootstrapping
#' nreps   = 496
#' vec.mu  = rep(0,nreps)
#' vec.var = rep(0,nreps)
#' for (i in 1:nreps){
#'    sample.id = boot.mblock(1000, b=50)
#'    sample.y  = vec.y[sample.id]
#'    vec.mu[i]  = mean(sample.y)
#'    vec.var[i] = var(sample.y)
#'    print(paste("iteration ",i,"/",nreps," complete.", sep=""))
#' }
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(vec.x, vec.y, type="l", main="1d signal")  # 1d signal
#' hist(vec.mu, main="mean CI", xlab="mu")         # mean
#' abline(v=data.mu, col="red", lwd=4)
#' hist(vec.var, main="variance CI", xlab="sigma") # variance
#' abline(v=data.var, col="blue", lwd=4)
#' par(opar)
#' 
#' @references 
#' \insertRef{kunsch_jackknife_1989}{maotai}
#' 
#' @export
boot.mblock <- function(N, b=max(2,round(N/10))){
  ###################################################################
  # Preprocessing
  myn = round(N)
  myb = round(b)
  vec1N = c(1:myn,1:myn,1:myn)
  
  ###################################################################
  # Preparation
  id0  = 1
  idb  = (myn-myb+1)
  id0b = (id0:idb) # starting point 
  
  ###################################################################
  # Computation
  output = c()
  while (length(output)<myn){
    idst   = base::sample(id0b,1)
    output = c(output,vec1N[idst:(idst+b-1)])
  }
  
  ###################################################################
  # Return
  return(output[1:myn])
}