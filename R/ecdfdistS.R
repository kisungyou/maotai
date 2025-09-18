#' Distance Measures between Samples through Empirical Cumulative Distribution Functions
#' 
#' 
#' We measure distance between two empirical cumulative distribution functions of the data. 
#' Unlike \code{\link[maotai]{ecdfdist}}, this function takes raw data samples as input, and 
#' internally computes the empirical cumulative distribution functions (ECDF) for distance calculations.
#' 
#' @param veclist a length \eqn{N} list of vectors.
#' @param method name of the distance/dissimilarity measure. Case insensitive (default: \code{ks}).
#' @param p exponent for \code{Lp} or \code{Wasserstein} distance (default: \code{p=1}).
#' @param as.dist a logical; \code{TRUE} to return \code{dist} object, \code{FALSE} to return an \eqn{(N\times N)} symmetric matrix of pairwise distances (default: \code{FALSE}).
#' 
#' @return either \code{dist} object of an \eqn{(N\times N)} symmetric matrix of pairwise distances by \code{as.dist} argument.
#' 
#' @examples
#' \donttest{
#' ## toy example : 10 of random and uniform distributions
#' mylist = list()
#' for (i in 1:10){
#'   mylist[[i]] = stats::rnorm(50, sd=2)
#' }
#' for (i in 11:20){
#'   mylist[[i]] = stats::runif(50, min=-5)
#' }
#' 
#' ## compute three distances
#' d_KS = ecdfdistS(mylist, method="KS")
#' d_LP = ecdfdistS(mylist, method="Lp")
#' d_OT = ecdfdistS(mylist, method="Wasserstein")
#' 
#' ## visualize
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(d_KS[,nrow(d_KS):1], axes=FALSE, main="Kolmogorov-Smirnov")
#' image(d_LP[,nrow(d_LP):1], axes=FALSE, main="Lp (p=1)")
#' image(d_OT[,nrow(d_OT):1], axes=FALSE, main="Wasserstein (p=1)")
#' par(opar)
#' }
#' 
#' @export
ecdfdistS <- function(veclist, method=c("KS","Lp","Wasserstein"), p=1, as.dist=FALSE){
  ###############################################
  # Preprocessing
  if (!is.list(veclist)){
    stop("* ecdfdistS : input 'veclist' should be a list of numeric vectors.")
  }
  if (!check_veclist(veclist)){
    stop("* ecdfdistS : input 'veclist' should be a list of numeric vectors.")
  }
  methodss = c("ks","wasserstein","lp")
  mymethod = tolower(method)
  mymethod = match.arg(mymethod, methodss)
  myp      = round(p)
  if (myp < 1){
    stop("* ecdfdistS : exponent 'p' should be >= 1.")
  }
  
  ################################################
  # Computation
  output = switch(mymethod,
                  "lp" = lp_ecdf_distance_matrix(veclist, myp),
                  "ks" = lp_ecdf_distance_matrix(veclist, Inf),
                  "wasserstein" = wasserstein_p_distance_matrix(veclist, myp))
  
  ###############################################
  # Report
  if (as.dist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}


# ## toy example : N(0,1), Unif[0,10], and N(-5,4)
# n_test  = 900
# n_third = round(n_test/3)
# 
# list_sam  = vector("list", length=n_test)
# list_ecdf = vector("list", length=n_test)
# 
# idnow = 0
# for (i in seq_len(n_third)){
#   idnow = idnow+1
#   n_sam = sample(100:500, 1)
#   list_sam[[idnow]] = stats::rnorm(n_sam, mean=0, sd=1)
#   list_ecdf[[idnow]] = stats::ecdf(list_sam[[idnow]])
# }
# for (i in seq_len(n_third)){
#   idnow = idnow+1
#   n_sam = sample(100:500, 1)
#   list_sam[[idnow]] = stats::runif(n_sam, min=0, max=10)
#   list_ecdf[[idnow]] = stats::ecdf(list_sam[[idnow]])
# }
# for (i in seq_len(n_third)){
#   idnow = idnow+1
#   n_sam = sample(100:500, 1)
#   list_sam[[idnow]] = stats::rnorm(n_sam, mean=-5, sd=2)
#   list_ecdf[[idnow]] = stats::ecdf(list_sam[[idnow]])
# }
# 
# R_ks = ecdfdist(list_ecdf, method="KS")
# R_lp = ecdfdist(list_ecdf, method="Lp")
# R_wt = ecdfdist(list_ecdf, method="Wasserstein")
# S_ks = ecdfdistS(list_sam, method="KS")
# S_lp = ecdfdistS(list_sam, method="Lp")
# S_wt = ecdfdistS(list_sam, method="Wasserstein")
# 
# par(mfrow=c(1,3), pty="s")
# image(abs(R_ks-S_ks), main=paste0("Diff:KS-",round(norm(R_ks-S_ks,"F"),3)))
# image(abs(R_lp-S_lp), main=paste0("Diff:LP-",round(norm(R_lp-S_lp,"F"),3)))
# image(abs(R_wt-S_wt), main=paste0("Diff:OT-",round(norm(R_wt-S_wt,"F"),3)))
# #   
# # 
# # 
# # 
# # # 
# # compute Kolmogorov-Smirnov distance
# microbenchmark::microbenchmark(
#   R_ks = ecdfdist(list_ecdf, method="KS"),
#   R_lp = ecdfdist(list_ecdf, method="Lp"),
#   R_wt = ecdfdist(list_ecdf, method="Wasserstein"),
#   S_ks = ecdfdistS(list_sam, method="KS"),
#   S_lp = ecdfdistS(list_sam, method="Lp"),
#   S_wt = ecdfdistS(list_sam, method="Wasserstein"),
#   times=5
# )
# 
# 
# 
# 
# 
