#' DP-means Algorithm for Clustering Euclidean Data
#' 
#' DP-means is a nonparametric clustering method motivated by DP mixture model in that 
#' the number of clusters is determined by a parameter \eqn{\lambda}. The larger 
#' the \eqn{\lambda} value is, the smaller the number of clusters is attained. 
#' In addition to the original paper, we added an option to randomly permute 
#' an order of updating for each observation's membership as a common 
#' heuristic in the literature of cluster analysis. 
#' 
#' @param data an \eqn{(n\times p)} data matrix for each row being an observation.
#' @param lambda a threshold to define a new cluster.
#' @param maxiter maximum number of iterations.
#' @param abstol stopping criterion 
#' @param permute.order a logical; \code{TRUE} if random order for permutation is used, \code{FALSE} otherwise.
#' 
#' @return a named list containing
#' \describe{
#' \item{cluster}{an \eqn{(n\times ndim)} matrix whose rows are embedded observations.}
#' \item{centers}{a list containing information for out-of-sample prediction.}
#' }
#' 
#' @examples
#' ## define data matrix of two clusters
#' x1  = matrix(rnorm(50*3,mean= 2), ncol=3)
#' x2  = matrix(rnorm(50*3,mean=-2), ncol=3)
#' X   = rbind(x1,x2)
#' lab = c(rep(1,50),rep(2,50))
#' 
#' ## run dpmeans with several lambda values
#' solA <- dpmeans(X, lambda= 5)$cluster
#' solB <- dpmeans(X, lambda=10)$cluster
#' solC <- dpmeans(X, lambda=20)$cluster
#' 
#' ## visualize the results
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,4), pty="s")
#' plot(X,col=lab,  pch=19, cex=.8, main="True", xlab="x", ylab="y")
#' plot(X,col=solA, pch=19, cex=.8, main="dpmeans lbd=5", xlab="x", ylab="y")
#' plot(X,col=solB, pch=19, cex=.8, main="dpmeans lbd=10", xlab="x", ylab="y")
#' plot(X,col=solC, pch=19, cex=.8, main="dpmeans lbd=20", xlab="x", ylab="y")
#' par(opar)
#' 
#' \donttest{
#' ## let's find variations by permuting orders of update
#' ## used setting : lambda=20, we will 8 runs
#' sol8 <- list()
#' for (i in 1:8){
#'   sol8[[i]] = dpmeans(X, lambda=20, permute.order=TRUE)$cluster
#' }
#' 
#' ## let's visualize
#' vpar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,4), pty="s")
#' for (i in 1:8){
#'   pm = paste("permute no.",i,sep="")
#'   plot(X,col=sol8[[i]], pch=19, cex=.8, main=pm, xlab="x", ylab="y")
#' }
#' par(vpar)
#' }
#' 
#' @references 
#' \insertRef{kulis_revisiting_2012}{maotai}
#' 
#' @export
dpmeans <- function(data, lambda=1, maxiter=1234, abstol=1e-6, permute.order=FALSE){
  ############################################################
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* dpmeans : an input 'data' should be a matrix without any missing/infinite values.")
  }
  # Parameter and Initialization
  n = nrow(data)
  p = ncol(data)
  k = 1                                   # set k=1
  labels = rep(1,n)                       # labels={1,2,...,n}
  mu     = matrix(colMeans(data), nrow=1) # global mean
  lambda = as.double(lambda)
  
  ############################################################
  # Main Iteration
  ss.old = compute.ss(data, labels, mu)
  ss.new = 0
  for (iter in 1:maxiter){
    # 0. updating order of observations
    if (permute.order){
      idseq = sample(1:n)
    } else {
      idseq = 1:n
    }
    # 1. update the class membership per each class
    for (i in idseq){
      # 1-1. compute distances to the centers
      # dic = rep(0, k); for (j in 1:k){dic[j] = sum((as.vector(data[i,])-as.vector(mu[j,]))^2)}
      dic = as.vector(dat2centers(data[i,], mu)); # cpp conversion
      # 1-2. assign new or stay
      if (min(dic) > lambda){
        k = k+1
        labels[i] = k
        mu = rbind(mu, data[i,])
      } else {
        idmins = which(dic==min(dic))
        if (length(idmins)>1){
          labels[i] = sample(idmins, 1)
        } else {
          labels[i] = idmins 
        }
      }
    }
    
    # 2. rearrange the label (remove empty ones)
    labels = as.factor(labels)
    ulabel = sort(unique(labels))
    labnew = rep(0,n)
    for (i in 1:length(ulabel)){
      labnew[(labels==ulabel[i])] = i
    }
    labels = labnew
    k      = round(max(labels))
    
    # 3. compute per-class means
    uassign = sort(unique(labels))
    mu = array(0,c(k,p))
    for (i in 1:k){
      idmean = which(labels==uassign[i])
      if (length(idmean)==1){
        mu[i,] = as.vector(data[idmean,])
      } else {
        mu[i,] = as.vector(colMeans(data[idmean,]))
      }
    }
    
    # 4. compute updated sum of squared errors
    ss.new   = compute.ss(data, labels, mu)
    ss.delta = ss.old-ss.new
    ss.old   = ss.new
    
    # 5. stop if updating is not significant
    if (ss.delta < abstol){
      break
    }
  }
  
  ############################################################
  # Return the results
  output = list()
  output$cluster = as.factor(labels)
  output$centers = mu
  return(output)
}

# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
compute.ss <- function(data, label, centers){
  p = ncol(data)
  if (is.vector(centers)){
    centers = matrix(centers, nrow=1)
  }
  ulabel = sort(unique(label))
  output = 0
  for (i in 1:length(ulabel)){
    subdata = data[(label==ulabel[i]),]
    if (!is.vector(subdata)){
      nn = nrow(subdata)
      for (j in 1:nn){
        output = output + sum((as.vector(subdata[j,])-as.vector(centers[i,]))^2)
      }
    }
  }
  return(output)
}
