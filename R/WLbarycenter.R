#' Barycenter of vMF Distributions Under a Wasserstein-Like Geometry
#' 
#' Given a collection of von Mises-Fisher (vMF) distributions, each characterized 
#' by a mean direction \eqn{\mathbf{\mu}} and a concentration parameter \eqn{\kappa},
#' this function solves the geometric mean problem to compute the barycentric vMF 
#' distribution under an approximate Wasserstein geometry.
#' 
#' @param means An \eqn{(n \times p)} matrix where each row represents the mean 
#' direction of one of the \eqn{n} vMF distributions.
#' @param concentrations A length-\eqn{n} vector of nonnegative concentration parameters.
#' @param weights A weight vector of length \eqn{n}. If \code{NULL}, equal weights 
#' (\code{rep(1/n, n)}) are used.
#' 
#' @return A named list containing:
#' \describe{
#'   \item{mean}{A length-\eqn{p} vector representing the barycenter direction.}
#'   \item{concentration}{A scalar representing the barycenter concentration.}
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
#' # Compute the barycenter under the Wasserstein-like geometry
#' barycenter <- WLbarycenter(means, concentrations)
#' 
#' # Convert barycenter mean direction to an angle
#' bary_angle <- atan2(barycenter$mean[2], barycenter$mean[1])
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' 
#' # Plot the unit circle
#' plot(cos(seq(0, 2 * pi, length.out = 200)), sin(seq(0, 2 * pi, length.out = 200)), 
#'      type = "l", col = "gray", lwd = 2, xlab = "x", ylab = "y", 
#'      main = "Barycenter of vMF Distributions on S^1")
#' 
#' # Add input mean directions
#' points(means[,1], means[,2], col = "blue", pch = 19, cex = 1.5)
#' 
#' # Add the computed barycenter
#' points(cos(bary_angle), sin(bary_angle), col = "red", pch = 17, cex = 2)
#' 
#' # Add legend
#' legend("bottomleft", legend = c("vMF Means", "Barycenter"), col = c("blue", "red"), 
#'        pch = c(19, 17), cex = 1)
#' 
#' # Plot the concentration parameters
#' hist(concentrations, main = "Concentration Parameters", xlab = "Concentration")
#' abline(v=barycenter$concentration, col="red", lwd=2)
#' par(opar)
#' }
#' 
#' @export
WLbarycenter <- function(means, concentrations, weights=NULL){
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
    stop("* WLbarycenter : cardinalities of the means and concentrations do not match.")
  }
  
  # weights
  if ((length(weights)==0)&&(is.null(weights))){
    data_weights = rep(1/length(data_kappa), length(data_kappa))
  } else {
    data_weights = as.vector(weights)
    data_weights = data_weights/base::sum(data_weights)
  }
  if (any(data_weights < .Machine$double.eps)||(length(data_weights)!=length(data_kappa))){
    stop("* WLbarycenter : invalid 'weights'. Please see the documentation.")
  }
  
  ###############################################
  # Computation
  # mean direction
  output_mean = as.vector(cpp_WL_weighted_mean(data_X, data_weights))
  
  # concentration
  kappa_tmp    = base::sum(data_weights/sqrt(data_kappa))
  kappa_output = 1/(kappa_tmp*kappa_tmp) 
  
  ###############################################
  # Return
  output = list(mean=output_mean,
                concentration=kappa_output)
  return(output)
}
