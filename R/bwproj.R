#' Bures-Wasserstein Projection
#' 
#' Projects a set of high-dimensional Gaussian distributions specified by their means and covariance matrices 
#' onto a lower-dimensional subspace using the distance-preserving projection method. 
#' 
#' @param means an \eqn{(n,p)} matrix of Gaussian means, where \eqn{n} is the number of Gaussians and \eqn{p} is the original dimension.
#' @param covs a \eqn{(p,p,n)} array of Gaussian covariance matrices.
#' @param target_dim an integer specifying the target lower dimension \eqn{d} (default is 2).
#' @param max_iter an integer specifying the maximum number of iterations for the optimization (default is 100).
#' @param verbose a logical flag indicating whether to print progress messages (default is TRUE).
#' 
#' @return a named list containing \describe{
#' \item{U}{the \eqn{(p,d)} projection matrix mapping original space to the lower-dimensional space.}
#' \item{proj_means}{the \eqn{(n,d)} matrix of projected Gaussian means.}
#' \item{proj_covs}{the \eqn{(d,d,n)} array of projected Gaussian covariance matrices.}
#' \item{iter_obj}{a vector of objective function values at each iteration.}
#' \item{iter_gnorm}{a vector of gradient norms at each iteration.}
#' }
#' 
#' @export
bwproj <- function(means, covs, target_dim=2, max_iter = 100, verbose=TRUE){
  # check
  if (!check_datamat(means)){
    stop("* bwproj: 'means' must be a matrix without NA/Inf.")
  }
  par_target_dim = max(1, round(target_dim))
  par_max_iter   = max(2, round(max_iter))
  par_verbose = as.logical(verbose)
  
  # run
  result = bwp_project_gaussians_rcpp(
    means, 
    covs, 
    d=par_target_dim, 
    max_iter=par_max_iter,
    verbose=par_verbose
  )
  
  # prepare
  d = par_target_dim
  p = dim(covs)[1]
  n = dim(covs)[3]
  
  output = list()
  output$U = result$U
  output$proj_means = means%*%output$U
  output$proj_covs = array(0,c(d,d,n))
  for (i in 1:n){
    output$proj_covs[,,i] = t(output$U)%*%covs[,,i]%*%output$U
  }
  output$iter_obj = result$obj
  output$iter_gnorm = result$grad_norm
  return(output)
}


# # personal example
# set.seed(123)
# 
# # Problem sizes
# K <- 3          # number of classes
# n_per <- 30     # Gaussians per class
# n <- K * n_per
# p <- 20
# d <- 2
# 
# # Class labels
# y <- rep(1:K, each = n_per)
# 
# # Means: n x p
# M <- matrix(0, n, p)
# 
# # Covariances: array dim=c(p,p,n)
# Sigma <- array(0, dim = c(p, p, n))
# 
# # Build three "prototype" covariance structures (SPD) + mean centers
# # Each class uses a different dominant subspace / conditioning pattern.
# make_spd_from_eigs <- function(eigs) {
#   Q <- qr.Q(qr(matrix(rnorm(p*p), p, p)))
#   S <- Q %*% diag(eigs) %*% t(Q)
#   (S + t(S))/2
# }
# 
# # Distinct eigenvalue profiles (controls shape/anisotropy)
# eigs1 <- c(rep(6, 3), rep(1, p-3))     # class 1: strong 3-d anisotropy
# eigs2 <- c(rep(1, 8), rep(4, 2), rep(1, p-10))  # class 2: different anisotropy
# eigs3 <- c(rep(2, p))                  # class 3: near-isotropic (but shifted mean)
# 
# S0_1 <- make_spd_from_eigs(eigs1)
# S0_2 <- make_spd_from_eigs(eigs2)
# S0_3 <- make_spd_from_eigs(eigs3)
# 
# # Mean centers (separated mostly in first few coordinates)
# sep <- 4.0   # try 3, 4, 5 ...
# 
# c1 <- sep * c(rep( 2.5, 3), rep(0, p-3))
# c2 <- sep * c(rep(-2.5, 3), rep(0, p-3))
# c3 <- sep * c(rep( 0.0, 3), rep(2.5, 3), rep(0, p-6))
# 
# # Within-class variability controls
# mean_noise <- 0.6
# cov_noise  <- 0.25   # magnitude of random SPD perturbation
# jitter     <- 0.15   # diagonal jitter for stability
# 
# for (i in 1:n) {
#   cls <- y[i]
# 
#   # Means: class center + noise
#   if (cls == 1) M[i,] <- c1 + rnorm(p, sd = mean_noise)
#   if (cls == 2) M[i,] <- c2 + rnorm(p, sd = mean_noise)
#   if (cls == 3) M[i,] <- c3 + rnorm(p, sd = mean_noise)
# 
#   # Covariances: class prototype + SPD perturbation
#   A <- matrix(rnorm(p*p), p, p)
#   P <- crossprod(A) / p  # SPD perturbation
#   if (cls == 1) S <- S0_1 + cov_noise * P
#   if (cls == 2) S <- S0_2 + cov_noise * P
#   if (cls == 3) S <- S0_3 + cov_noise * P
# 
#   Sigma[,,i] <- (S + t(S))/2 + jitter * diag(p)
# }
# 
# # run
# # Run projection
# res <- bwproj(
#   means = M,
#   covs = Sigma,
#   target_dim = d,
#   max_iter = 200,
# )
# 
# Uhat <- res$U
# 
# # Draw a 2D Gaussian confidence ellipse (base graphics)
# # mu: length-2 mean
# # S : 2x2 SPD covariance
# # level: confidence level (0.95 by default)
# draw_gaussian_ellipse <- function(mu, S, level = 0.95, npoints = 200,
#                                   border = "black", lwd = 1, lty = 1) {
#   stopifnot(length(mu) == 2, all(dim(S) == c(2,2)))
#   S <- (S + t(S))/2
# 
#   # Chi-square radius for 2D
#   r2 <- qchisq(level, df = 2)
#   r  <- sqrt(r2)
# 
#   # eigendecomp for ellipse
#   ee <- eigen(S, symmetric = TRUE)
#   vals <- pmax(ee$values, 0)
#   vecs <- ee$vectors
# 
#   theta <- seq(0, 2*pi, length.out = npoints)
#   unit  <- rbind(cos(theta), sin(theta))  # 2 x npoints
# 
#   # Transform unit circle -> ellipse
#   A <- vecs %*% diag(sqrt(vals), 2, 2)
#   pts <- sweep(A %*% unit, 1, mu, "+")    # 2 x npoints
# 
#   lines(pts[1,], pts[2,], col = border, lwd = lwd, lty = lty)
#   invisible(pts)
# }
# 
# plot_projected_gaussians <- function(M, Sigma, Uhat, y, level = 0.95,
#                                      show_centers = TRUE,
#                                      cols = c("tomato", "royalblue", "darkgreen"),
#                                      lwd = 1.5, lty = 1) {
#   stopifnot(is.matrix(M), length(dim(Sigma)) == 3)
#   n <- nrow(M)
#   p <- ncol(M)
#   stopifnot(dim(Sigma)[1] == p, dim(Sigma)[2] == p, dim(Sigma)[3] == n)
#   stopifnot(nrow(Uhat) == p, ncol(Uhat) == 2)
#   stopifnot(length(y) == n)
# 
#   # Embedded means and covariances
#   M_low <- M %*% Uhat  # n x 2
#   Sig_low <- array(0, dim = c(2,2,n))
#   for (i in 1:n) {
#     Sig_low[,,i] <- t(Uhat) %*% Sigma[,,i] %*% Uhat
#     Sig_low[,,i] <- (Sig_low[,,i] + t(Sig_low[,,i]))/2
#   }
# 
#   # Determine plot limits from ellipses (cheap bounding)
#   # Use the largest axis length for each ellipse: sqrt(qchisq)*sqrt(max eigenvalue)
#   r <- sqrt(qchisq(level, df = 2))
#   rad <- numeric(n)
#   for (i in 1:n) {
#     ev <- eigen(Sig_low[,,i], symmetric = TRUE, only.values = TRUE)$values
#     rad[i] <- r * sqrt(max(pmax(ev, 0)))
#   }
#   xlim <- range(M_low[,1] + c(-1,1) %o% rad)
#   ylim <- range(M_low[,2] + c(-1,1) %o% rad)
# 
#   # Base plot canvas
#   plot(NA, xlim = xlim, ylim = ylim, asp = 1,
#        xlab = "Projected coordinate 1", ylab = "Projected coordinate 2",
#        main = sprintf("Projected Gaussians (%d%% ellipses)", round(100*level)))
# 
#   # Draw ellipses
#   for (i in 1:n) {
#     cls <- y[i]
#     border_col <- cols[cls]
#     draw_gaussian_ellipse(mu = M_low[i,], S = Sig_low[,,i],
#                           level = level, border = border_col, lwd = lwd, lty = lty)
#     if (show_centers) {
#       points(M_low[i,1], M_low[i,2], pch = 19, cex = 0.6, col = border_col)
#     }
#   }
# 
#   legend("topright", legend = paste("Class", sort(unique(y))),
#          col = cols[sort(unique(y))], lwd = lwd, bty = "n")
# 
#   invisible(list(M_low = M_low, Sigma_low = Sig_low))
# }
# 
# # After running:
# # res <- bwproj(...)
# # Uhat <- res$U
# 
# out <- plot_projected_gaussians(
#   M = M,
#   Sigma = Sigma,
#   Uhat = Uhat,
#   y = y,
#   level = 0.95,
#   show_centers = TRUE
# )
