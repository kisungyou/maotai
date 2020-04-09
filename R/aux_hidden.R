# Hidden Functions for Future Use
# these functions can be loaded using 'utils::getFromNamespace'
# by the command 'getFromNamespace("function_name","maotai");
# 
# Here, all the functions require 'diss' object from 'stats' package.
#
# 01. hidden_kmedoids : PAM algorithm 
# 02. hidden_bmds     : Bayesian Multidimensional Scaling
# 03. hidden_cmds     : Classical Multidimensional Scaling
# 04. hidden_kmeanspp : k-means++ algorithm.
# 05. hidden_tsne     : t-SNE visualization.




# 01. hidden_kmedoids -----------------------------------------------------
#' @keywords internal
hidden_kmedoids <- function(distobj, nclust=2){
  myk = round(nclust)
  return(cluster::pam(distobj, k = myk))
}

# 02. hidden_bmds ---------------------------------------------------------
#' @keywords internal
hidden_bmds <- function(x, ndim=2, par.a=5, par.alpha=0.5, par.step=1, mc.iter=8128, verbose=FALSE){
  ######################################################
  # Initialization
  ndim   = round(ndim)
  embedy = hidden_cmds(x, ndim)$embed
  x      = as.matrix(x)
  ndim   = round(ndim)
  
  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* bmds - 'ndim' should be an integer in [1,ncol(data)). ")
  }
  
  n = nrow(x)
  m = n*(n-1)/2
  
  ######################################################
  # Preliminary Computation
  # 1. apply CMDS for initialization
  y     = as.matrix(base::scale(embedy,       # (N x ndim) centered 
                                center=TRUE, scale=FALSE)) 
  Delta = as.matrix(stats::dist(y))           # (N x N) pairwise distances
  
  # 2. initialization
  eigy   = base::eigen(stats::cov(y))       
  X0     = y%*%eigy$vectors         # (N x ndim) rotated
  gamma0 = diag(X0)                 # variances ?
  sigg0  = compute_SSR(x, Delta)/m;  
  beta0  = apply(X0,2,var)/2
  
  # 3. run the main part
  runcpp <- main_bmds(x, X0, sigg0, par.a, par.alpha, mc.iter, par.step, verbose, beta0)
  Xsol   <- runcpp$solX
  Xdist  <- as.matrix(stats::dist(Xsol))
  
  output = list()
  output$embed  = Xsol
  output$stress = compute_stress(x, Xdist)
  return(output)
  # return Rcpp::List::create(Rcpp::Named("solX")=Xsol,Rcpp::Named("solSSR")=SSRsol);
}


# 03. hidden_cmds ---------------------------------------------------------
#' @keywords internal
hidden_cmds <- function(x, ndim=2){
  ##################################################3
  # Check Input and Transform
  ndim = round(ndim)
  k  = as.integer(ndim)
  x  = as.matrix(x)
  D2 = (x^2)          # now squared matrix
  n  = nrow(D2)
  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* cmds - 'ndim' should be an integer in [1,ncol(data)). ")
  }
  
  ##################################################3
  # Computation
  J = diag(n) - (1/n)*outer(rep(1,n),rep(1,n))
  B = -0.5*J%*%D2%*%J
  eigB = eigen(B)
  
  LL = eigB$values[1:k]
  EE = eigB$vectors[,1:k]
  
  # Y  = as.matrix(base::scale((EE%*%diag(sqrt(LL))), center=TRUE, scale=FALSE))
  Y  = EE%*%diag(sqrt(LL))
  DY = as.matrix(stats::dist(Y))
  
  output = list()
  output$embed  = Y
  output$stress = compute_stress(x, DY)
  return(output)
}

# 04. hidden_kmeanspp -----------------------------------------------------
#' @keywords internal
hidden_kmeanspp <- function(x, k=2){
  ##################################################3
  # Check Input and Transform
  x  = as.matrix(x)
  n  = nrow(x)
  K  = round(k)
  if (K >= n){
    stop("* kmeanspp : 'k' should be smaller than the number of observations.")
  }
  if (K < 2){
    stop("* kmeanspp : 'k' should be larger than 1.")
  }
  id.now = 1:n
  
  ##################################################3
  # Computation
  #   initialize
  id.center = base::sample(id.now, 1)
  id.now    = base::setdiff(id.now, id.center)
  #   iterate
  for (i in 1:(K-1)){
    # compute distance to the nearest
    tmpdmat = x[id.now, id.center]
    if (i==1){
      d2vec = as.vector(tmpdmat)^2
      d2vec = d2vec/base::sum(d2vec)
    } else {
      d2vec = as.vector(base::apply(tmpdmat, 1, base::min))^2
      d2vec = d2vec/base::sum(d2vec)
    }
    # sample one
    id.tmp = base::sample(id.now, 1, prob=d2vec)
    # update
    id.center = c(id.center, id.tmp)
    id.now    = base::setdiff(id.now, id.tmp)
  }
  #   let's compute label
  dmat    = x[,id.center]
  cluster = base::apply(dmat, 1, base::which.min)
  
  ##################################################
  # Return
  return(cluster)
}


# 05. hidden_tsne ---------------------------------------------------------
#' @keywords internal
hidden_tsne <- function(dx, ndim=2, ...){
  ##################################################
  # Pass to 'Rtsne'
  k      = round(ndim)
  dx     = as.matrix(dx)
  tmpout = Rtsne::Rtsne(dx, dims=k, ..., is_distance=TRUE)
  Y   = tmpout$Y
  DY  = as.matrix(stats::dist(Y))
  
  ##################################################
  # Return
  output = list()
  output$embed  = Y
  output$stress = compute_stress(dx, DY)
  return(output)
}