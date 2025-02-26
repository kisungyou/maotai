#' Extract meaningful information from the von Mises-Fisher mixture model
#' 
#' Given a mixture of von Mises-Fisher distributions, this function computes 
#' several related quantities of the data on the unit hypersphere with respect 
#' to the specified model.
#' 
#' @param data an \eqn{(n\times d)} data matrix.
#' @param means an \eqn{(k\times d)} matrix of means. 
#' @param concentrations a vector of length \eqn{k} of concentration parameters.
#' @param weights a vector of length \eqn{k} of mixing weights.
#' 
#' @return a named list containing \describe{
#' \item{densities}{a vector of length \eqn{n} of the densities of the data points.}
#' \item{clustering}{a vector of length \eqn{n} of the hard clustering results.}
#' \item{loglkd}{the log-likelihood of the data.}
#' \item{AIC}{the Akaike information criterion.}
#' \item{BIC}{the Bayesian information criterion.}
#' }
#' 
#' @export
movMF_info <- function(data, means, concentrations, weights){
  ###############################################
  # Preprocessing
  # data
  if (!is.matrix(data)){
    my_data = cpp_WL_normalise(as.matrix(data))
  } else {
    my_data = cpp_WL_normalise(data)
  }
  
  # means
  if (!is.matrix(means)){
    my_means = cpp_WL_normalise(as.matrix(means))
  } else {
    my_means = cpp_WL_normalise(means)
  }
  
  # concentrations
  my_concentrations = as.vector(concentrations)
  if (length(my_concentrations)!=base::nrow(my_means)){
    stop("* movMF_info : cardinalities of the means and concentrations do not match.")
  }
  
  # weights
  if ((length(weights)==0)&&(is.null(weights))){
    my_weights = rep(1/length(my_concentrations), length(my_concentrations))
  } else {
    my_weights = as.vector(weights)
    my_weights = my_weights/base::sum(my_weights)
  }
  if (any(my_weights < .Machine$double.eps)||(length(my_weights)!=length(my_concentrations))){
    stop("* movMF_info : invalid 'weights'. Please see the documentation.")
  }
  
  ###############################################
  # BASIC COMPUTATION
  par_n = base::nrow(my_data)
  par_p = base::ncol(my_data)
  par_k = base::length(my_weights)
  membership <- array(0,c(par_n, par_k))
  for (k in 1:par_k){
    membership[,k] <- aux_vmf_density(my_data, my_means[k,], my_concentrations[k])*my_weights[k]
  }
  
  ###############################################
  # ADVANCED
  # density (n,) vector
  out_densities = base::rowSums(membership) 
  
  # hard clustering results
  out_cluster = rep(0, par_n)
  for (i in 1:par_n){
    tmp_clust <- as.vector(membership[i,])
    tmp_clust <- tmp_clust/base::sum(tmp_clust)
    out_cluster[i] <- which.max(tmp_clust)
  }
  
  # log-likelihood
  out_loglkd <- base::sum(base::log(out_densities))
  
  # information criteria
  df <- (par_k - 1) + par_k*((par_p-1)+1)
  out_aic <- -2*out_loglkd + 2*df
  out_bic <- -2*out_loglkd + base::log(par_n)*df
  
  ###############################################
  # RETURN
  output = list()
  output$densities = out_densities
  output$clustering = out_cluster
  output$loglkd <- out_loglkd
  output$AIC <- out_aic 
  output$BIC <- out_bic
  
  return(output)
}



# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
aux_vmf_density <- function(x, mu, kappa){
  # x : (n x d) data matrix
  # mu : (d,)   mean vector
  # kappa : >0  concentration parameter
  
  # initialization
  d <- base::ncol(x)
  
  # normalizing constant
  log_c_p_kappa <- (d/2 - 1) * log(kappa) - (d/2) * log(2 * pi) - log(base::besselI(kappa, nu = d/2 - 1, expon.scaled = TRUE)) - kappa
  c_p_kappa <- exp(log_c_p_kappa)
  
  
  # iterate
  densities = rep(0, base::nrow(x))
  for (i in 1:base::nrow(x)){
    densities[i] <- base::exp(kappa*base::sum(as.vector(x[i,])*mu))*c_p_kappa
  }
  return(densities)
  
}

# # simple example
# # data matrix normalized
# X = as.matrix(iris[,2:4])
# X = as.matrix(scale(X, center=TRUE, scale=FALSE))
# X = X/sqrt(rowSums(X^2))
# 
# # fit the model with movMF package
# out_movMF <- movMF::movMF(X, 3)
# clust_movMF <- predict(out_movMF, X)
# 
# # use my function
# fit_weights <- out_movMF$alpha
# fit_means   <- out_movMF$theta/sqrt(rowSums(out_movMF$theta^2))
# fit_concentrations <- sqrt(rowSums(out_movMF$theta^2))
# clust_mine <- movMF_info(X, fit_means, fit_concentrations, fit_weights)$clustering
