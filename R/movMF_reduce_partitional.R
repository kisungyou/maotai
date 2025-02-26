#' von Mises-Fisher mixture model reduction - partitional method
#' 
#' When given parameters of the von Mises-Fisher mixture model, this function 
#' aims at mixture model reduction using a partitional method.
#' 
#' @param means a \eqn{(K \times p)} matrix of means of the von Mises-Fisher components.
#' @param concentrations a \eqn{K} vector of concentrations of the von Mises-Fisher components.
#' @param weights a \eqn{K} vector of weights of the von Mises-Fisher components.
#' @param target.num a desired number of components after reduction. Default is 2.
#' @param method a clustering method to be used. Default is "hclust".
#' 
#' @return a named list of the reduced mixture model containing \describe{
#' \item{means}{a \eqn{(\code{target.num} \times p)} matrix of means of the von Mises-Fisher components.}
#' \item{concentrations}{a \eqn{\code{target.num}} vector of concentrations of the von Mises-Fisher components.}
#' \item{weights}{a \eqn{\code{target.num}} vector of weights of the von Mises-Fisher components.}
#' }
#' 
#' @export
movMF_reduce_partitional <- function(means, concentrations, weights, target.num=2, method=c("hclust","kmedoids")){
  ###############################################
  # Preprocessing
  # means
  if (!is.matrix(means)){
    data_means = cpp_WL_normalise(as.matrix(means))
  } else {
    data_means = cpp_WL_normalise(means)
  }
  d = base::ncol(data_means)
  
  # concentrations
  data_concentrations = as.vector(concentrations)
  if (length(data_concentrations)!=base::nrow(data_means)){
    stop("* reduce_movMF_greedy : cardinalities of the means and concentrations do not match.")
  }
  
  # weights
  if ((length(weights)==0)&&(is.null(weights))){
    data_weights = rep(1/length(data_concentrations), length(data_concentrations))
  } else {
    data_weights = as.vector(weights)
    data_weights = data_weights/base::sum(data_weights)
  }
  if (any(data_weights < .Machine$double.eps)||(length(data_weights)!=length(data_concentrations))){
    stop("* movMF_reduce_partitional : invalid 'weights'. Please see the documentation.")
  }
  
  # desired number
  K = length(data_weights)
  par_target_num = round(target.num)
  if ((par_target_num <= 1)||(par_target_num >= K)){
    stop("* movMF_reduce_partitional : 'target.num' must be greater than 1 and less than the number of components.")
  }
  
  ###############################################
  # MAIN ROUTINE
  # check the clustering method
  par_clustering = match.arg(method)
  
  # pairwise distance computation
  pdist_mat = array(0,c(K,K))
  for (i in 1:(K-1)){
    for (j in (i+1):K){
      term1 = sphere_dist(as.vector(data_means[i,]), as.vector(data_means[j,]))
      term2sq = (d-1)*((1/sqrt(data_concentrations[i]) - 1/sqrt(data_concentrations[j]))^2)
      pdist_mat[i,j] <- pdist_mat[j,i] <- sqrt((term1*term1) + term2sq)
    }
  }
  pdist_obj <- stats::as.dist(pdist_mat)
  
  # clustering
  if (par_clustering=="hclust"){
    clust_obj <- fastcluster::hclust(pdist_obj, method="single")
    obtained_clust <- stats::cutree(clust_obj, k=par_target_num)
  } else {
    clust_obj <- cluster::pam(pdist_obj, k=par_target_num)
    obtained_clust <- clust_obj$clustering
  }
  
  # get the indices
  list_indices = vector("list", length=par_target_num)
  for (i in 1:par_target_num){
    list_indices[[i]] = which(obtained_clust==i)
  }
  
  # for each indices, compute the barycenter
  new_means = array(0,c(par_target_num, base::nrow(data_means)))
  new_concentrations = rep(0, par_target_num)
  new_weights = rep(0, par_target_num)
  
  for (it in 1:par_target_num){
    # current indices
    now_index = list_indices[[it]]
    
    # branching : single vs multiple observations in a cluster
    if (length(now_index)==1){
      new_means[it,] = data_means[now_index,]
      new_concentrations[it] = data_concentrations[now_index]
      new_weights[it] = data_weights[now_index]
    } else {
      tmp_means = data_means[now_index,]
      tmp_concentration = data_concentrations[now_index]
      tmp_weights = data_weights[now_index]
      new_weights[it] = base::sum(tmp_weights) # aggregated weights
      tmp_weights = tmp_weights/base::sum(tmp_weights)
      
      # compute the barycenter
      tmp_bary <- WLbarycenter(tmp_means, tmp_concentration, tmp_weights)
      new_means[it,] <- tmp_bary$mean
      new_concentrations[it] <- tmp_bary$concentration
    }
  }
  
  # return the output
  output = list(means=new_means,
                concentrations=new_concentrations,
                weights=new_weights)
  return(output)
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
# red_hclust <- movMF_reduce_partitional(big_means, big_concentrations, big_weights, target.num=3, method="hclust")
# red_medoid <- movMF_reduce_partitional(big_means, big_concentrations, big_weights, target.num=3, method="kmedoids")
# 
# clust_hclust <- movMF_info(X, red_hclust$means, red_hclust$concentrations, red_hclust$weights)$clustering
# clust_medoid <- movMF_info(X, red_medoid$means, red_medoid$concentrations, red_medoid$weights)$clustering
# 
# # visualize
# par(mfrow=c(1,3), pty="s")
# plot(X, col=clust_big_movMF, pch=19, main="Original")
# plot(X, col=clust_hclust, pch=19, main="Reduced (hclust)")
# plot(X, col=clust_medoid, pch=19, main="Reduced (kmedoids)")
