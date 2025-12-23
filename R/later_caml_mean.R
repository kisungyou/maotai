#' CAML: Frechet mean 
#' 
#' 
#' 
#' 
#' 
#' @keywords internal
#' @noRd
caml_mean <- function(list_samples, maxiter=100, abstol=1e-6){
  # check the input
  if (!check_datalist(list_samples)){
    stop("* caml_mean: 'list_samples' must have the same-sized matrices.")
  }
  par_maxiter = max(5, round(maxiter))
  par_abstol = max(1e-12, as.double(abstol))
  
  # weights
  N = length(list_samples)
  w = rep(1/N, N)
  
  # compute
  cpp_run <- caml_frechet_mean(
    Y_list_R = list_samples, # uncentered version is fine
    w = w,
    maxit = par_maxiter, 
    tol      = par_abstol,
    step0    = 1.0,
    c1       = 1e-4,
    beta     = 0.5,
    verbose  = TRUE
  )
 
  return(cpp_run)
}

# set.seed(42)
# 
# 
# # "True" latent positions (unknown in real scenarios)
# 
# mlbench_true = mlbench::mlbench.smiley()
# Y_true  <- as.matrix(scale(mlbench_true$x, center=TRUE, scale=FALSE))
# Y_class <- mlbench_true$classes
# 
# n <- nrow(Y_true)
# k <- ncol(Y_true)
# S <- 200  # number of MCMC samples
# 
# 
# # Generate MCMC-like samples: rotated + noisy
# list_samples <- vector("list", S)
# 
# for (s in 1:S) {
#   Qs <- qr.Q(qr(matrix(rnorm(k*k), k, k)))  # random orthogonal
#   noise <- 0.1 * matrix(rnorm(n*k), n, k)
#   list_samples[[s]] <- Y_true %*% Qs + noise
# }
# 
# # Compute the CAML Fréchet mean
# res <- caml_mean(list_samples, maxiter = 200, abstol = 1e-8)
# 
# Y_hat  <- res$Y
# X_hat  <- res$YYT
# value  <- res$cost
# 
# # plot
# par(mfrow=c(1,3), pty="s")
# plot(Y_true, col=Y_class, main="True shape")
# plot(Y_hat, col=Y_class, main="Frechet mean")
# image(X_hat, main="Mean Gram matrix")
# 
# 
# # caml stands for
# # Coherent Analysis of MCMC samples from LPM models
# # present this as a framework
