#' EP-means Algorithm for Clustering Empirical Distributions
#' 
#' EP-means is a variant of k-means algorithm adapted to cluster 
#' multiple empirical cumulative distribution functions under metric structure 
#' induced by Earth Mover's Distance.
#' 
#' @param elist a length \eqn{N} list of either vector or \code{ecdf} objects.
#' @param k the number of clusters.
#' 
#' @return a named list containing \describe{
#' \item{cluster}{an integer vector indicating the cluster to which each \code{ecdf} is allocated.}
#' \item{centers}{a length \eqn{k} list of centroid \code{ecdf} objects.}
#' }
#' 
#' @examples
#' \donttest{
#' ## two sets of 1d samples, 10 each and add some noise
#' #    set 1 : mixture of two gaussians
#' #    set 2 : single gamma distribution
#' 
#' # generate data
#' elist = list()
#' for (i in 1:10){
#'    elist[[i]] = stats::ecdf(c(rnorm(100, mean=-2), rnorm(50, mean=2)))
#' }
#' for (j in 11:20){
#'    elist[[j]] = stats::ecdf(rgamma(100,1) + rnorm(100, sd=sqrt(0.5)))
#' }
#' 
#' # run EP-means with k clusters 
#' # change the value below to see different settings
#' myk   = 2
#' epout = epmeans(elist, k=myk)
#' 
#' # visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,myk))
#' for (k in 1:myk){
#'   idk = which(epout$cluster==k)
#'   for (i in 1:length(idk)){
#'     if (i<2){
#'       pm = paste("class ",k," (size=",length(idk),")",sep="")
#'       plot(elist[[idk[i]]], verticals=TRUE, lwd=0.25, do.points=FALSE, main=pm)
#'     } else {
#'       plot(elist[[idk[i]]], add=TRUE, verticals=TRUE, lwd=0.25, do.points=FALSE)
#'     }
#'     plot(epout$centers[[k]], add=TRUE, verticals=TRUE, lwd=2, col="red", do.points=FALSE)
#'   }
#' }
#' par(opar)
#' }
#' 
#' @references 
#' \insertRef{henderson_epmeans_2015}{maotai}
#' 
#' @export
epmeans <- function(elist, k=2){
  ###############################################
  # Preprocessing
  clist = elist_epmeans(elist)  # will use quantized ones only / flist = elist_fform(qlist)
  myk   = round(k)
  myn   = length(clist)
  
  # Quantization
  mylength = 1000
  qseq = seq(from=1e-6, to=1-(1e-6), length.out=mylength)
  qmat = array(0,c(myn,mylength))
  for (n in 1:myn){
    qmat[n,] = as.vector(stats::quantile(clist[[n]], qseq))
  }
  
  ###############################################
  # Rcpp k-means
  tmpcpp = cpp_kmeans(qmat, myk)$means
  
  ###############################################
  # Pairwise Distance Computation
  #   wrap
  mylist1 = list()
  mylist2 = list()
  for (n in 1:myn){
    mylist1[[n]] = stats::ecdf(as.vector(qmat[n,]))
  }
  for (k in 1:myk){
    mylist2[[k]] = stats::ecdf(as.vector(tmpcpp[k,]))
  }
  #   compute pairwise distance using Earth Mover's Distance
  pdistmat = dist2_wasserstein(mylist1, mylist2, 1)
  #   index
  label = base::apply(pdistmat, 1, which.min)
  
  
  ###############################################
  # Return : we want to add 'Silhouette'
  output = list()
  output$cluster = as.integer(label)
  output$centers = mylist2
  return(output)
}
  

# ## personal examples
# cdf0  = stats::ecdf(rnorm(100, sd=3))    # original ECDF
# qseq  = seq(from=0,to=1,length.out=1000) # quantile sequence
# quant = stats::quantile(cdf0, qseq)
# cdf1  = stats::ecdf(quant)
# 
# par(mfrow=c(1,2))
# plot(cdf0, main="Original")
# plot(cdf1, main="Recovered")


