#' Geometric Median of vMF Distributions Under a Wasserstein-Like Geometry
#' 
#' Given a collection of von Mises-Fisher (vMF) distributions, each characterized 
#' by a mean direction \eqn{\mathbf{\mu}} and a concentration parameter \eqn{\kappa},
#' this function solves the geometric median problem to compute the vMF 
#' distribution that minimizes the weighted sum of distances under an approximate Wasserstein geometry.
#' 
#' @param means An \eqn{(n \times p)} matrix where each row represents the mean 
#' direction of one of the \eqn{n} vMF distributions.
#' @param concentrations A length-\eqn{n} vector of nonnegative concentration parameters.
#' @param weights A weight vector of length \eqn{n}. If \code{NULL}, equal weights 
#' (\code{rep(1/n, n)}) are used.
#' 
#' @return A named list containing:
#' \describe{
#'   \item{mean}{A length-\eqn{p} vector representing the median direction.}
#'   \item{concentration}{A scalar representing the median concentration.}
#' }
#' 
#' @examples
#' \donttest{
#' # Set seed for reproducibility
#' set.seed(123)
#' 
#' # Number of vMF distributions
#' n <- 5   
#' 
#' # Generate mean directions concentrated around a specific angle (e.g., 45 degrees)
#' base_angle <- pi / 4  # 45 degrees in radians
#' angles <- rnorm(n, mean = base_angle, sd = pi / 20)  # Small deviation from base_angle
#' means <- cbind(cos(angles), sin(angles))  # Convert angles to unit vectors
#' 
#' # Generate concentration parameters with large magnitudes (tight distributions)
#' concentrations <- rnorm(n, mean = 50, sd = 5)  # Large values around 50
#' 
#' # Compute the median under the Wasserstein-like geometry
#' barycenter <- WLmedian(means, concentrations)
#' 
#' # Convert median mean direction to an angle
#' bary_angle <- atan2(barycenter$mean[2], barycenter$mean[1])
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' 
#' # Plot the unit circle
#' plot(cos(seq(0, 2 * pi, length.out = 200)), sin(seq(0, 2 * pi, length.out = 200)), 
#'      type = "l", col = "gray", lwd = 2, xlab = "x", ylab = "y", 
#'      main = "Median of vMF Distributions on S^1")
#' 
#' # Add input mean directions
#' points(means[,1], means[,2], col = "blue", pch = 19, cex = 1.5)
#' 
#' # Add the computed barycenter
#' points(cos(bary_angle), sin(bary_angle), col = "red", pch = 17, cex = 2)
#' 
#' # Add legend
#' legend("bottomleft", legend = c("vMF Means", "Median"), col = c("blue", "red"), 
#'        pch = c(19, 17), cex = 1)
#' 
#' # Plot the concentration parameters
#' hist(concentrations, main = "Concentration Parameters", xlab = "Concentration")
#' abline(v=barycenter$concentration, col="red", lwd=2)
#' par(opar)
#' }
#' 
#' 
#' @export
WLmedian <- function(means, concentrations, weights=NULL){
  ###############################################
  # Preprocessing
  # means
  if (!is.matrix(means)){
    data_X = cpp_WL_normalise(as.matrix(means))
  } else {
    data_X = cpp_WL_normalise(means)
  }
  
  # concentrations
  data_kappa = as.vector(concentrations)
  if (length(data_kappa)!=base::nrow(data_X)){
    stop("* WLmedian : cardinalities of the means and concentrations do not match.")
  }
  
  # weights
  if ((length(weights)==0)&&(is.null(weights))){
    data_weights = rep(1/length(data_kappa), length(data_kappa))
  } else {
    data_weights = as.vector(weights)
    data_weights = data_weights/base::sum(data_weights)
  }
  if (any(data_weights < .Machine$double.eps)||(length(data_weights)!=length(data_kappa))){
    stop("* WLmedian : invalid 'weights'. Please see the documentation.")
  }
  
  
  ###############################################
  # COMPUTATION
  # initialize
  old_mean  = as.vector(cpp_WL_weighted_mean(data_X, data_weights))
  old_kappa = base::mean(data_kappa) 
  
  # iterate
  par_nsample = base::nrow(data_X)
  par_dim     = base::ncol(data_X)
  for (it in 1:100){
    # compute the new weights
    now_weights = rep(0, par_nsample)
    for (i in 1:par_nsample){
      delta <- old_mean - as.vector(data_X[i, ])
      c <- sqrt(sum(delta^2))                # chord length
      term1 <- (2 * asin(pmin(1, c/2)))^2    # squared geodesic
      term2 = (par_dim)*((1/sqrt(old_kappa) - 1/sqrt(data_kappa[i]))^2)
      dist_old2sam   = sqrt(term1+term2)
      now_weights[i] = data_weights[i]/dist_old2sam 
    }
    now_weights = now_weights/base::sum(now_weights)
    
    # compute the new mean and kappa
    new_mean  = as.vector(cpp_WL_weighted_mean(data_X, now_weights))
    new_w_sum = base::sum(now_weights)
    new_kappa_tmp = base::sum(now_weights/sqrt(data_kappa))
    new_kappa = (new_w_sum/new_kappa_tmp)^2
    
    # update
    inc_mean  = sqrt(sum((old_mean-new_mean)^2))
    inc_kappa = abs(old_kappa-new_kappa)
    inc_max   = max(inc_mean, inc_kappa)
    
    old_mean = new_mean
    old_kappa = new_kappa
    if (inc_max < 1e-8){
      break
    }
  }
  
  ###############################################
  # Return
  output = list(mean=old_mean,
                concentration=old_kappa)
  return(output)
}
