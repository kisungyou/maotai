# Hidden Functions for Future Use
# these functions can be loaded using 'utils::getFromNamespace'
# by the command 'getFromNamespace("function_name","maotai");
# 
# Here, all the functions require 'diss' object from 'stats' package.
#
# 01. hidden_kmedoids      : PAM algorithm 
#     hidden_kmedoids_best : PAM algorithm + use Silhouette (maximum)
# 02. hidden_bmds          : Bayesian Multidimensional Scaling
# 03. hidden_cmds          : Classical Multidimensional Scaling
# 04. hidden_kmeanspp      : k-means++ algorithm.
# 05. hidden_tsne          : t-SNE visualization.
# 06. hidden_nem           : Negative Eigenvalue Magnitude
# 07. hidden_nef           : Negative Eigenfraction
# 08. hidden_emds          : Euclified Multidimensional Scaling
# 09. hidden_hclust        : FASTCLUSTER - hclust function
# 10. hidden_dbscan        : DBSCAN      - dbscan function
# 11. hidden_silhouette    : mimics that of cluster's silhouette
# 12. hidden_mmds          : metric multidimensional scaling by SMACOF 
# 13. hidden_PHATE         : return row-stochastic matrix & time stamp



# 00. hidden_checker ------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_checker <- function(xobj){
  if (inherits(xobj, "dist")){
    return(as.matrix(xobj))
  } else if (inherits(xobj, "matrix")){
    check1 = (nrow(xobj)==ncol(xobj))
    check2 = isSymmetric(xobj)
    check3 = all(abs(diag(xobj))<.Machine$double.eps*10)
    if (check1&&check2&&check3){
      return(as.matrix(xobj))
    } else {
      stop("* hidden : matrix is not valid.")
    }
  } else {
    stop("* hidden : input is not valid.")
  }
}

# 01. hidden_kmedoids & hidden_kmedoids_best ------------------------------
#' @keywords internal
#' @noRd
hidden_kmedoids <- function(distobj, nclust=2){
  distobj = stats::as.dist(hidden_checker(distobj))
  myk     = round(nclust)
  return(cluster::pam(distobj, k = myk))
}
#' @keywords internal
#' @noRd
hidden_kmedoids_best <- function(distobj, mink=2, maxk=10){
  # prepare
  kvec = seq(from=round(mink),to=round(maxk), by = 1)
  knum = length(kvec)
  svec = rep(0,knum)
  nobj = nrow(as.matrix(distobj))
  clut = array(0,c(nobj,knum))
  
  for (k in 1:knum){
    know = kvec[k]
    if (know < 2){
      svec[k]  = 0
      clut[,k] = rep(1,nobj)
    } else {
      pamx     = hidden_kmedoids(distobj, nclust=kvec[k])
      svec[k]  = pamx$silinfo$avg.width
      clut[,k] = pamx$clustering  
    }
  }
  # return the output
  output = list()
  output$opt.k = kvec[which.max(svec)]
  output$score = svec     # knum-vector of silhouette scores
  output$label = clut     # (n,knum) cluster labels
  return(output)
}

# 02. hidden_bmds ---------------------------------------------------------
#' @keywords internal
#' @noRd
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
#' @noRd
hidden_cmds <- function(x, ndim=2){
  ##################################################3
  # Check Input and Transform
  x    = hidden_checker(x)
  ndim = round(ndim)
  k  = as.integer(ndim)
  D2 = (x^2)          # now squared matrix
  n  = nrow(D2)
  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* cmds : 'ndim' should be an integer in [1,ncol(data)). ")
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
#' @noRd
hidden_kmeanspp <- function(x, k=2){
  ##################################################3
  # Check Input and Transform
  x  = hidden_checker(x)
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
  output = list()
  output$center  = id.center
  output$cluster = cluster
  return(output)
}


# 05. hidden_tsne ---------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_tsne <- function(dx, ndim=2, ...){
  ##################################################
  # Pass to 'Rtsne'
  dx     = hidden_checker(dx)
  k      = round(ndim)
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

# 06. hidden_nem ----------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_nem <- function(xdiss){
  ##################################################3
  # Check Input and Transform
  xx = hidden_checker(xdiss)
  D2 = (xx^2)
  n  = nrow(D2)

  ##################################################3
  # Computation
  H = diag(n) - (1/n)*base::outer(rep(1,n),rep(1,n))
  S = -0.5*(H%*%D2%*%H)
  eigS = base::eigen(S, only.values = TRUE)
  evals = eigS$values

  ##################################################3
  # Finalize
  output = abs(min(evals))/max(evals)
  return(output)
}


# 07. hidden_nef ----------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_nef <- function(xdiss){
  ##################################################3
  # Check Input and Transform
  xx = hidden_checker(xdiss)
  D2 = (xx^2)
  n  = nrow(D2)
  
  ##################################################3
  # Computation
  H = diag(n) - (1/n)*base::outer(rep(1,n),rep(1,n))
  S = -0.5*(H%*%D2%*%H)
  eigS = base::eigen(S)
  evals = eigS$values
  
  ##################################################3
  # Finalize
  output = sum(abs(evals[which(evals<0)]))/sum(abs(evals))
  return(output)
}

# 08. hidden_emds ---------------------------------------------------------
#' Euclified Multidimensional Scaling
#' 
#' strategy 1 : transitive closure of the triangle inequality (labdsv)
#' strategy 2 : Non-Euclidean or Non-metric Measures Can Be Informative; adding positive numbers to all off-diagonal entries
#' 
#' @keywords internal
#' @noRd
hidden_emds <- function(xdiss, ndim=2, method=c("closure","gram")){
  ##################################################3
  # Check Input and Transform
  x    = hidden_checker(xdiss)
  ndim = round(ndim)
  k  = as.integer(ndim)
  n  = nrow(x)
  if ((length(ndim)>1)||(ndim<1)||(ndim>=nrow(x))){
    stop("* emds : 'ndim' should be an integer in [1,nrow(x)). ")
  }
  
  method = match.arg(method) 
  mydim  = round(ndim)
  
  ##################################################3
  # Branching
  if (hidden_nef(x) < 100*.Machine$double.eps){ # if Euclidean, okay
    output = hidden_cmds(x, ndim=mydim)
  } else { # if not Euclidean
    if (method=="closure"){ # strategy 1 : transitive closure of the triangle inequality
      xnew   = as.matrix(labdsv::euclidify(stats::as.dist(x))) # well it seems to work well..
      output = hidden_cmds(xnew, ndim = mydim)
    } else {                # strategy 2 : add positive numbers to all off-diagonal entries
      gamma0 = emds_gamma0(x)
      ggrid  = seq(from=min(0.001, gamma0/1000), to=(gamma0*0.999), length.out=20) # just try 20 cases
      vgrid  = rep(0,20)
      for (i in 1:20){
        xtmp = x + ggrid[i]
        diag(xtmp) = rep(0,nrow(xtmp))
        vgrid[i]   = hidden_nef(xtmp)
      }
      idopts = which.min(vgrid)
      if (length(idopts)>1){ # if multiple, use the first one.
        idopts = idopts[1]
      }
      optgamma   = ggrid[idopts]
      xnew       = x + optgamma
      diag(xnew) = rep(0,nrow(xnew))
      output     = hidden_cmds(xnew, ndim = mydim)
    } 
  }
  
  ##################################################3
  # Report 
  return(output)
}

# 09. hidden_hclust -------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_hclust <- function(xdiss, mymethod, mymembers){
  return(fastcluster::hclust(xdiss, method=mymethod,
                             members=mymembers))
}

# 10. hidden_dbscan -------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_dbscan <- function(Xdiss, myeps, myminPts=5, ...){
  output = dbscan::dbscan(Xdiss, eps = myeps, minPts=myminPts, ...)
  return(output)
}

# 11. hidden_silhouette --------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_silhouette <- function(xdiss, label){
  x    = as.integer(as.factor(label))
  hsil = cluster::silhouette(x, xdiss)
  
  output = list()
  output$local  = as.vector(hsil[,3])
  output$global = base::mean(as.vector(hsil[,3]))
  return(output)
}

# 12. hidden_mmds          : metric multidimensional scaling by SMACOF  --------
#' @keywords internal
#' @noRd
hidden_mmds <- function(x, ndim=2, maxiter=200, abstol=1e-5){
  # Check Input and Transform
  x    = hidden_checker(x)
  ndim = round(ndim)
  myiter = max(50, round(maxiter))
  mytol  = max(100*.Machine$double.eps, as.double(abstol))
  
  # Run with Rcpp
  return(cpp_mmds(x, ndim, myiter, mytol))
}



# 13. hidden_PHATE --------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_PHATE <- function(x, nbdk=5, alpha=2){
  # Check Input and Transform
  x = hidden_checker(x) # now it's a matrix
  n = base::nrow(x)
  nbdk  = max(1, round(nbdk))
  alpha = max(sqrt(.Machine$double.eps), as.double(alpha))
  
  # k-th nearest distance
  nndist = rep(0,n)
  for (i in 1:n){
    tgt = as.vector(x[i,])
    nndist[i] = tgt[order(tgt)[nbdk+1]]
  }
  
  # Build Kernel Matrix
  matK = array(1,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      term1 = exp(-((x[i,j]/nndist[i])^(alpha)))
      term2 = exp(-((x[i,j]/nndist[j])^(alpha)))
      matK[i,j] <- matK[j,i] <- 0.5*(term1+term2)
    }
  }
  vecD = base::rowSums(matK)
  matP = base::diag(1/vecD)%*%matK
  matA = base::diag(1/base::sqrt(vecD))%*%matK%*%base::diag(1/base::sqrt(vecD))
  
  # Eigenvalues and Von-Neumann Entropy
  eigA  = eigen(matA)$values
  eigA  = eigA[(eigA>0)]
  vec.t = 1:1000
  vec.H = rep(0,1000)
  for (i in 1:1000){
    eig.t  = eigA^i
    eig.t  = eig.t/base::sum(eig.t)
    term.t = -base::sum(eig.t*base::log(eig.t))
    if (is.na(term.t)){
      vec.t = vec.t[1:(i-1)]
      vec.H = vec.H[1:(i-1)]
      break
    } else {
      vec.H[i] = term.t
    }
  }
  
  # Optimal Stopping Criterion
  Pout  = matP
  opt.t = round(hidden_knee_clamped(vec.t, vec.H))
  for (i in 1:(opt.t-1)){
    Pout = Pout%*%matP
  }
  
  # return the output
  output = list()
  output$P = Pout
  output$t = opt.t
  return(output)
}

# X   = as.matrix(iris[,1:4])
# lab = as.factor(iris[,5])
# 
# D    = stats::dist(X)
# cmd2 = cmdscale(D, k=2)
# mmdA = hidden_mmds(D, ndim=2, abstol=1e-2)
# mmdB = hidden_mmds(D, ndim=2, abstol=1e-10)
# 
# par(mfrow=c(1,3), pty="s")
# plot(cmd2, col=lab, main = "cmds")
# plot(mmdA, col=lab, main="mmds-2")
# plot(mmdB, col=lab, main="mmds-8")


# # example -----------------------------------------------------------------
# library(labdsv)
# data(bryceveg) # returns a vegetation data.frame
# dis.bc <- as.matrix(dsvdis(bryceveg,'bray/curtis')) # calculate a Bray/Curtis
# 
# emds = getFromNamespace("hidden_emds","maotai")
# out.cmds <- cmds(dis.bc, ndim=2)$embed
# out.emds1 <- emds(dis.bc, ndim=2, method="closure")$embed
# out.emds2 <- emds(dis.bc, ndim=2, method="gram")$embed
# 
# par(mfrow=c(3,1))
# plot(out.cmds, main="cmds")
# plot(out.emds1, main="emds::closure")
# plot(out.emds2, main="emds::gram")



