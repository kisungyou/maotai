#' Large-scale Generalized Procrustes Analysis
#' 
#' We modify generalized Procrustes analysis for large-scale data by 
#' first setting a subset of anchor points and applying the attained transformation 
#' to the rest data. If \code{sub.id} is a vector \code{1:dim(x)[1]}, it uses all 
#' observations as anchor points, reducing to the conventional generalized Procrustes analysis.
#' 
#' @param x a \eqn{(k\times m\times n)} 3d array, where \eqn{k} is the number of points, \eqn{m} the number of dimensions, and \eqn{n} the number of samples.
#' @param sub.id a vector of indices for defining anchor points.
#' @param scale a logical; \code{TRUE} if scaling is applied, \code{FALSE} otherwise.
#' 
#' @return a \eqn{(k\times m\times n)} 3d array of aligned samples.
#' 
#' @examples 
#' \dontrun{
#' ## This should be run if you have 'shapes' package installed.
#' library(shapes)
#' data(gorf.dat)
#' 
#' ## apply anchor-based method and original procGPA
#' out.proc = shapes::procGPA(gorf.dat, scale=TRUE)$rotated # procGPA from shapes package
#' out.anc4 = lgpa(gorf.dat, sub.id=c(1,4,9,7), scale=TRUE) # use 4 points 
#' out.anc7 = lgpa(gorf.dat, sub.id=1:7, scale=TRUE)        # use all but 1 point as anchors
#' 
#' ## visualize
#' opar = par(mfrow=c(3,4), pty="s")
#' plot(out.proc[,,1], main="procGPA No.1", pch=18)
#' plot(out.proc[,,2], main="procGPA No.2", pch=18)
#' plot(out.proc[,,3], main="procGPA No.3", pch=18)
#' plot(out.proc[,,4], main="procGPA No.4", pch=18)
#' plot(out.anc4[,,1], main="4 Anchors No.1", pch=18, col="blue")
#' plot(out.anc4[,,2], main="4 Anchors No.2", pch=18, col="blue")
#' plot(out.anc4[,,3], main="4 Anchors No.3", pch=18, col="blue")
#' plot(out.anc4[,,4], main="4 Anchors No.4", pch=18, col="blue")
#' plot(out.anc7[,,1], main="7 Anchors No.1", pch=18, col="red")
#' plot(out.anc7[,,2], main="7 Anchors No.2", pch=18, col="red")
#' plot(out.anc7[,,3], main="7 Anchors No.3", pch=18, col="red")
#' plot(out.anc7[,,4], main="7 Anchors No.4", pch=18, col="red")
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{goodall_procrustes_1991}{maotai}
#' 
#' @author Kisung You
#' @export
lgpa <- function(x, sub.id = 1:(dim(x)[1]), scale=TRUE){
  ###################################################################
  # check : x
  if ((!is.array(x))||(length(dim(x))!=3)){
    stop("* lgpa : input 'x' should be a 3d array.")
  }
  dimsx = dim(x)
  k = dimsx[1]
  m = dimsx[2]
  n = dimsx[3]
  
  # check : sub.id
  sub.id = round(sub.id)
  sub.id = base::intersect(sub.id, 1:k)
  if ((max(sub.id) > k)||(!is.vector(sub.id))){
    stop("* lgpa : an input 'sub.id' should be a vector containing indices in [1,nrow(x)].")
  }
  
  ###################################################################
  # computation
  #   1. select out the subarray and compute means
  nsubid = length(sub.id)
  xsub   = x[sub.id,,]
  meanvecs = list()
  for (i in 1:n){
    meanvecs[[i]] = colMeans(xsub[,,i])
  }
  for (i in 1:n){
    xsub[,,i] = xsub[,,i] - matrix(rep(meanvecs[[i]],nsubid), ncol=m, byrow = TRUE)
  }
  #   2. compute PGA
  if (scale){
    xout = shapes::procGPA(xsub, scale=TRUE)$rotated  
  } else {
    xout = shapes::procGPA(xsub, scale=FALSE)$rotated  
  }
  
  #   3. compute rotation matrix
  rotmats = list()
  for (i in 1:n){
    tgt1 = xsub[,,i]
    tgt2 = xout[,,i]
    rotmats[[i]] = base::solve(t(tgt1)%*%tgt1, t(tgt1)%*%tgt2)
  }
  #   4. final computation
  output = array(0,dim(x))
  for (i in 1:n){
    tgtx = x[,,i]
    output[,,i] = (tgtx - matrix(rep(meanvecs[[i]],k), ncol=m, byrow = TRUE))%*%(rotmats[[i]])
  }
  
  ###################################################################
  # Report
  return(output)
}
