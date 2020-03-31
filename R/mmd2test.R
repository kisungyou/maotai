#' Kernel Two-sample Test with Maximum Mean Discrepancy
#' 
#' Maximum Mean Discrepancy (MMD) as a measure of discrepancy between 
#' samples is employed as a test statistic for two-sample hypothesis test 
#' of equal distributions. Kernel matrix \eqn{K} is a symmetric square matrix 
#' that is positive semidefinite.
#' 
#' @param K kernel matrix or an object of \code{kernelMatrix} class from \pkg{kernlab} package.
#' @param label label vector of class indices. 
#' @param method type of estimator to be used. \code{"b"} for biased and \code{"u"} for unbiased estimator of MMD.
#' @param mc.iter the number of Monte Carlo resampling iterations.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided kernel matrix.}
#' }
#' 
#' @examples 
#' ## small test for CRAN submission
#' dat1 <- matrix(rnorm(60, mean= 1), ncol=2) # group 1 : 30 obs of mean  1
#' dat2 <- matrix(rnorm(50, mean=-1), ncol=2) # group 2 : 25 obs of mean -1
#' 
#' dmat <- as.matrix(dist(rbind(dat1, dat2)))  # Euclidean distance matrix
#' kmat <- exp(-(dmat^2))                      # build a gaussian kernel matrix
#' lab  <- c(rep(1,30), rep(2,25))             # corresponding label
#' 
#' mmd2test(kmat, lab)                         # run the code !
#' 
#' \dontrun{
#' ## WARNING: computationally heavy. 
#' #  Let's compute empirical Type 1 error at alpha=0.05
#' niter = 496  
#' pvals1 = rep(0,niter)
#' pvals2 = rep(0,niter)
#' for (i in 1:niter){
#'   dat = matrix(rnorm(200),ncol=2)
#'   lab = c(rep(1,50), rep(2,50))
#'   lbd = 0.1
#'   kmat = exp(-lbd*(as.matrix(dist(dat))^2))
#'   pvals1[i] = mmd2test(kmat, lab, method="b")$p.value
#'   pvals2[i] = mmd2test(kmat, lab, method="u")$p.value
#'   print(paste("iteration ",i," complete..",sep=""))
#' }
#' 
#' #  Visualize the above at multiple significance levels
#' alphas = seq(from=0.001, to=0.999, length.out=100)
#' errors1 = rep(0,100)
#' errors2 = rep(0,100)
#' for (i in 1:100){
#'    errors1[i] = sum(pvals1<=alphas[i])/niter
#'    errors2[i] = sum(pvals2<=alphas[i])/niter
#' }
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(alphas, errors1, "b", main="Biased Estimator Error", 
#'      xlab="alpha", ylab="error", cex=0.5)
#' abline(a=0,b=1, lwd=1.5, col="red")
#' plot(alphas, errors2, "b", main="Unbiased Estimator Error", 
#'      xlab="alpha", ylab="error", cex=0.5)
#' abline(a=0,b=1, lwd=1.5, col="blue")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{gretton_kernel_2012}{maotai}
#' 
#' @export
mmd2test <- function(K, label, method=c("b","u"), mc.iter=999){
  ###############################################
  # Preprocessing
  DNAME = deparse(substitute(K))
  # 1. K : kernel matrix
  if (inherits(K, "kernelMatrix")){
    kmat = as.matrix(K)
  } else {
    kmat = as.matrix(K)
  }
  cond1 = (is.matrix(kmat))
  cond2 = (nrow(K)==ncol(kmat))
  cond3 = isSymmetric(kmat)
  if (!(cond1&&cond2&&cond3)){
    stop("* mmd2test : 'K' should be a kernel matrix.")
  }
  mineval = min(base::eigen(kmat, only.values = TRUE)$values)
  if (mineval<0){
    wm = paste("* mmd2test : 'K' may not be PD. Minimum eigenvalue is ",mineval,".",sep="")
    warning(wm)
  }
  # 2. label
  label = as.vector(as.integer(as.factor(label)))
  if ((length(label)!=nrow(kmat))||(length(unique(label))!=2)){
    stop("* mmd2test : 'label' should be a vector of proper length with 2 classes.")
  }
  ulabel = unique(label)
  # 3. method
  allmmm = c("b","u")
  method = match.arg(tolower(method), allmmm)
  
  ###############################################
  # compute statistic
  id1 = which(label==ulabel[1]); m = length(id1)
  id2 = which(label==ulabel[2]); n = length(id2)
  thestat = switch(method,
                   "b" = mmd_biased(kmat[id1,id1], kmat[id2,id2], kmat[id1,id2]),
                   "u" = mmd_unbiased(kmat[id1,id1], kmat[id2,id2], kmat[id1,id2]))
  
  ###############################################
  # Iteration
  mciter = round(mc.iter)
  itervals = rep(0,mciter)
  for (i in 1:mciter){
    permuted = sample(m+n)
    tmpid1 = permuted[1:m]
    tmpid2 = permuted[(m+1):(m+n)]
    
    itervals[i] = switch(method,
                         "b" = mmd_biased(kmat[tmpid1,tmpid1], kmat[tmpid2,tmpid2], kmat[tmpid1,tmpid2]),
                         "u" = mmd_unbiased(kmat[tmpid1,tmpid1], kmat[tmpid2,tmpid2], kmat[tmpid1,tmpid2]))
  }
  pvalue = (sum(itervals>=thestat)+1)/(mciter+1)
  
  ###############################################
  # REPORT
  hname   = "Kernel Two-sample Test with Maximum Mean Discrepancy"
  Ha    = "two distributions are not equal"
  names(thestat) = "MMD"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}



# compute two squared statistics ------------------------------------------
#' @keywords internal
#' @noRd
mmd_biased <- function(XX, YY, XY){
  # parameters
  m = nrow(XX)
  n = nrow(YY)
  
  # computation
  return((sum(XX)/(m^2)) + (sum(YY)/(n^2)) - ((2/(m*n))*sum(XY)))
}
#' @keywords internal
#' @noRd
mmd_unbiased <- function(XX, YY, XY){
  # parameters
  m = nrow(XX)
  n = nrow(YY)
  
  # computation
  term1 = (sum(XX)-sum(diag(XX)))/(m*(m-1))
  term2 = (sum(YY)-sum(diag(YY)))/(n*(n-1))
  term3 = (2/(m*n))*sum(XY)
  return((term1+term2-term3))
}
