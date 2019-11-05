# auxiliary functions for computation -------------------------------------
# (1) aux_pinv       : pseudo-inverse
# (2) aux_pseudomean : compute distance from 1st observation to pseudo mean by rest points
 
# (1) aux_pinv ------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_pinv <- function(A){
  svdA      = base::svd(A)
  tolerance = (.Machine$double.eps)*max(c(nrow(A),ncol(A)))*as.double(max(svdA$d))
  
  idxcut    = which(svdA$d <= tolerance)
  invDvec   = (1/svdA$d)
  invDvec[idxcut] = 0
  
  output = (svdA$v%*%diag(invDvec)%*%t(svdA$u))
  return(output)
}


# (2) aux_pseudomean ------------------------------------------------------
#' @keywords internal
#' @noRd
aux_pseudomean <- function(dmat){
  # we need embedding .. umm .. automatic dimension selection
  if (nrow(dmat)==1){
    stop("* aux_pseudomean : error..")
  } else if (nrow(dmat)==2){
    return(dmat[1,2])
  } else {
    embedded = aux_pseudomean_auto(dmat)
    n = nrow(embedded)
    p = ncol(embedded)
    
    # centering based on other points
    emcenter = as.vector(base::colMeans(embedded[2:n,]))
    embednew = embedded - matrix(rep(emcenter,n), ncol=p, byrow=TRUE)
    
    # compute scalar
    d1mat = dmat[2:n,2:n]                          # d(x,y)
    d2mat = as.matrix(stats::dist(embednew[2:n,])) # ||x-y||
    d12mat = (d1mat*d2mat)
    d22mat = (d2mat^2)
    dlower = base::lower.tri(d12mat)
    cstar =sum(d12mat[dlower])/sum(d22mat[dlower])
    
    # update embednew and compute 
    erow1 = cstar*as.vector(embednew[1,])
    return(sqrt(sum(erow1^2))) 
  }
}
#' @keywords internal
#' @noRd
aux_pseudomean_auto <- function(dmat){ # only positive eigenvalues' part
  n = nrow(dmat)
  J = diag(rep(1,n))-(1/n)*outer(rep(1,n),rep(1,n))
  B = -(J%*%(dmat^2)%*%J)/2.0
  eigB = base::eigen(B, symmetric = TRUE) # decreasing order
  
  m = max(length(which(eigB$values > 0)),2)
  X = (eigB$vectors[,1:m])%*%(base::diag(sqrt(eigB$values[1:m])))
  return(X)
}

# # personal test : seems like it's working well enough
# x = rnorm(5, mean=3)
# y = matrix(rnorm(10*5),ncol=5)
# 
# dmat = as.matrix(dist(rbind(x,y)))
# val.alg = aux_pseudomean(dmat)
# val.true = sqrt(sum((x-as.vector(colMeans(y)))^2))
