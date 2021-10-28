# Hidden Functions : Computation
# 
# (01) hidden_gaussbary_2002R
#      - Ruschendorf, Uckelmann (2002) : On the n-coupling problem algorithm
#
# (02) hidden_gaussbary_2016A
#      - Alvarez-Esteban et al (2016) : A fixed-point approach to barycenters in Wasserstein space
#      - See Theorem 4.2 for the updating rules.
#
# (03) hidden_projsimplex_sorting : method by Wang and Carreira-Perpinan
#      Held, P. Wolfe, and H. Crowder, “Validation of subgradient optimization,” Mathematical Programming, vol. 6, pp. 62–88, 1974
#      (above is known by https://www.optimization-online.org/DB_FILE/2014/08/4498.pdf)
#      url : https://eng.ucmerced.edu/people/wwang5/papers/SimplexProj.pdf


# (01) hidden_gaussbary_2002R ---------------------------------------------
#      array3d : (p x p x N) array of SPD matrices.
#      weight  : N vector of weights. C++ handles L1 normalization.
#
#      RETURN a list containing
#             $mean : (p x p) barycenter matrix
#             $iter : the number of iterations
#' @keywords internal
#' @noRd
hidden_gaussbary_2002R <- function(array3d, weight, maxiter=50, abstol=(1e-8)){
  # minimal check
  weight = as.vector(weight)
  if (length(weight)!=dim(array3d)[3]){
    stop("* hidden_gaussbary_2002R : non-matching weight.")
  }
  
  # compute
  return(src_gaussbary_2002R(array3d, weight, maxiter, abstol))
}

# (02) hidden_gaussbary_2016A ---------------------------------------------
#      array3d : (p x p x N) array of SPD matrices.
#      weight  : N vector of weights. C++ handles L1 normalization.
#
#      RETURN a list containing
#             $mean : (p x p) barycenter matrix
#             $iter : the number of iterations
#' @keywords internal
#' @noRd
hidden_gaussbary_2016A <- function(array3d, weight, maxiter=50, abstol=(1e-8)){
  # minimal check
  weight = as.vector(weight)
  if (length(weight)!=dim(array3d)[3]){
    stop("* hidden_gaussbary_2016A : non-matching weight.")
  }
  
  # compute
  return(src_gaussbary_2016A(array3d, weight, maxiter, abstol))
}



# # test scenario 1 : random near-identity covariances
# # test scenario 2 : generate randomly some data + add noise
# p = 10
# n = 50
# noise = 0.01
# 
# somedat = matrix(runif(50*p), ncol=p)
# mycovs = array(0, c(p,p,n))
# for (i in 1:n){
#   mycovs[,,i] = cov(somedat + matrix(rnorm(50*p), ncol=p)*noise)
#   # print(paste0("The ",i,"-th matrix has rank=",round(Matrix::rankMatrix(mycovs[,,i]))))
# }
# myweight = rep(1/n, n)
# 
# 
# res_2002 = hidden_gaussbary_2002R(mycovs, myweight)
# res_2016 = hidden_gaussbary_2016A(mycovs, myweight)
# par(mfrow=c(1,3), pty="s")
# image(cov(somedat), main="true covariance")
# image(res_2002$mean, main=paste0("02R:iter=",round(res_2002$iter),"/error=",round(norm(res_2002$mean-cov(somedat),"F"),3)))
# image(res_2016$mean, main=paste0("16A:iter=",round(res_2016$iter),"/error=",round(norm(res_2016$mean-cov(somedat),"F"),3)))
# 
# microbenchmark(res_2002 = hidden_gaussbary_2002R(mycovs, myweight),
#                res_2016 = hidden_gaussbary_2016A(mycovs, myweight))



# (03) hidden_projsimplex_sorting ---------------------------------------------
#' @keywords internal
#' @noRd
hidden_projsimplex_sorting <- function(input){
  # preprocessing
  y = as.vector(input)
  D = length(y)
  u = base::sort(y, decreasing = TRUE)
  
  # main part
  uvec = rep(0,D)
  for (j in 1:D){
    uvec[j] = u[j] + (1/j)*(1-base::sum(u[1:j]))
  }
  rho = max(which(uvec > 0))
  lbd = (1/rho)*(1 - base::sum(u[1:rho]))
  
  # finalize
  x = pmax(y+lbd, 0)
  return(x)
}
# xx = matrix(rnorm(1000*2), ncol=2)
# yy = array(0,c(1000,2))
# for (i in 1:1000){
#   yy[i,] = hidden_projsimplex_13W(xx[i,])
# }
# plot(xx, main="scatter")
# points(yy, pch=19, col="red")
