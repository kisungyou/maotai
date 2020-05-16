#' Generate Index for Stationary Bootstrapping
#' 
#' Assuming data being dependent with cardinality \code{N}, \code{boot.stationary} returns 
#' a vector of index that is used for stationary bootstrapping. To describe, starting points 
#' are drawn from uniform distribution over \code{1:N} and the size of each block is 
#' determined from geometric distribution with parameter \eqn{p}.
#' 
#' 
#' @param N the number of observations.
#' @param p parameter for geometric distribution with the size of each block.
#' 
#' @return a vector of length \code{N} for moving block bootstrap sampling.
#' 
#' @examples
#' \donttest{
#' ## example : bootstrap confidence interval of mean and variances
#' vec.x = seq(from=0,to=10,length.out=100)
#' vec.y = sin(1.21*vec.x) + 2*cos(3.14*vec.x) + rnorm(100,sd=1.5)
#' data.mu  = mean(vec.y)
#' data.var = var(vec.y)
#' 
#' ## apply stationary bootstrapping
#' nreps   = 50
#' vec.mu  = rep(0,nreps)
#' vec.var = rep(0,nreps)
#' for (i in 1:nreps){
#'    sample.id = boot.stationary(100)
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
#' }
#' 
#' @references 
#' \insertRef{politis_stationary_1994}{maotai}
#' 
#' @export
boot.stationary <- function(N, p=0.25){
  ###################################################################
  # Preprocessing
  myn   = round(N)
  myp   = as.double(p)
  vec1n = 1:myn
  
  ###################################################################
  # Computation
  output = c()
  while (length(output)<myn){
    # 1. select the starting index
    idst = base::sample(vec1n,1)
    # 2. select the number
    nb = 0
    while (nb < 1){
      nb = stats::rgeom(1, myp)
    }
    # 3. select indices
    cnow = idst:(idst+nb-1)
    cnow = (cnow%%myn)
    cnow[(cnow==0)] = myn
    # 4. concatenate
    output = c(output, cnow)
  }
  
  ###################################################################
  # Return
  return(output[1:myn])
}