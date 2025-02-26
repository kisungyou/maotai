#' von Mises-Fisher mixture model reduction - Greedy method
#' 
#' When given parameters of the von Mises-Fisher mixture model, this function 
#' aims at mixture model reduction using a greedy method.
#' 
#' @param means a \eqn{(K \times p)} matrix of means of the von Mises-Fisher components.
#' @param concentrations a \eqn{K} vector of concentrations of the von Mises-Fisher components.
#' @param weights a \eqn{K} vector of weights of the von Mises-Fisher components.
#' @param target.num a desired number of components after reduction. Default is 2.
#' 
#' @return a named list of the reduced mixture model containing \describe{
#' \item{means}{a \eqn{(\code{target.num} \times p)} matrix of means of the von Mises-Fisher components.}
#' \item{concentrations}{a \eqn{\code{target.num}} vector of concentrations of the von Mises-Fisher components.}
#' \item{weights}{a \eqn{\code{target.num}} vector of weights of the von Mises-Fisher components.}
#' }
#' @export
movMF_reduce_greedy <- function(means, concentrations, weights, target.num=2){
  ###############################################
  # Preprocessing
  # means
  if (!is.matrix(means)){
    data_means = cpp_WL_normalise(as.matrix(means))
  } else {
    data_means = cpp_WL_normalise(means)
  }
  
  # concentrations
  data_concentrations = as.vector(concentrations)
  if (length(data_concentrations)!=base::nrow(data_means)){
    stop("* movMF_reduce_greedy : cardinalities of the means and concentrations do not match.")
  }
  
  # weights
  if ((length(weights)==0)&&(is.null(weights))){
    data_weights = rep(1/length(data_concentrations), length(data_concentrations))
  } else {
    data_weights = as.vector(weights)
    data_weights = data_weights/base::sum(data_weights)
  }
  if (any(data_weights < .Machine$double.eps)||(length(data_weights)!=length(data_concentrations))){
    stop("* movMF_reduce_greedy : invalid 'weights'. Please see the documentation.")
  }
  
  # desired number
  K = length(data_weights)
  par_target_num = round(target.num)
  if ((par_target_num <= 1)||(par_target_num >= K)){
    stop("* movMF_reduce_greedy : 'target.num' must be greater than 1 and less than the number of components.")
  }
  
  
  ###############################################
  # MAIN ROUTINE
  old_vmf <- list(means=data_means, concentrations=data_concentrations, weights=data_weights)
  
  while (length(old_vmf$weights) > par_target_num){
    # compute the pairwise distances
    dists = maotai::WLpdist(old_vmf$means, old_vmf$concentrations)
      
    # find the closest pair
    min_dist = Inf
    min_i = 0
    min_j = 0
    for (i in 1:(length(old_vmf$weights)-1)){
      for (j in (i+1):length(old_vmf$weights)){
        if (dists[i,j] < min_dist){
          min_dist = dists[i,j]
          min_i = i
          min_j = j
        }
      }
    }
    
    # save the non-closest pair
    new_means = old_vmf$means[-c(min_i, min_j),]
    new_concentrations = old_vmf$concentrations[-c(min_i, min_j)]
    new_weights = old_vmf$weights[-c(min_i, min_j)]
    
    # merge the closest pairs
    tgt_means = old_vmf$means[c(min_i, min_j),]
    tgt_kappa = old_vmf$concentrations[c(min_i, min_j)]
    tgt_alpha = old_vmf$weights[c(min_i, min_j)]
    tgt_alpha_norm = tgt_alpha/base::sum(tgt_alpha)
    tgt_barycenter <- maotai::WLbarycenter(tgt_means, tgt_kappa, tgt_alpha_norm)
    
    # update the old one
    old_vmf = list(means=rbind(new_means, tgt_barycenter$mean), 
                   concentrations=c(new_concentrations, tgt_barycenter$concentration),
                   weights=c(new_weights, base::sum(tgt_alpha)))
  }
  
  # return
  return(old_vmf)
}


# # simple example
# # data matrix normalized
# X = as.matrix(iris[,2:4])
# X = as.matrix(scale(X, center=TRUE, scale=FALSE))
# X = X%*%eigen(cov(X))$vectors[,1:2] # apply PCA
# X = X/sqrt(rowSums(X^2))
# 
# # fit the model with movMF package
# big_movMF <- movMF::movMF(X, 10)
# clust_big_movMF <- predict(big_movMF, X)
# convert_movMF <- movMF_convert(big_movMF)
# 
# big_means <- convert_movMF$means
# big_weights <- convert_movMF$weights
# big_concentrations <- convert_movMF$concentrations
# 
# # reduce to 3 components using different methods
# red3 <- movMF_reduce_greedy(big_means, big_concentrations, big_weights, target.num=3)
# red5 <- movMF_reduce_greedy(big_means, big_concentrations, big_weights, target.num=5)
# 
# clust3 <- movMF_info(X, red3$means, red3$concentrations, red3$weights)$clustering
# clust5 <- movMF_info(X, red5$means, red5$concentrations, red5$weights)$clustering
# 
# # visualize
# par(mfrow=c(1,3), pty="s")
# plot(X, col=clust_big_movMF, pch=19, main="Original")
# plot(X, col=clust3, pch=19, main="Greedy K=3")
# plot(X, col=clust5, pch=19, main="Greedy K=5")
