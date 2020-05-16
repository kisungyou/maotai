#' Trace Ratio Optimation
#' 
#' This function provides several algorithms to solve the following problem
#' \deqn{\textrm{max} \frac{tr(V^\top A V)}{tr(V^\top B V)} \textrm{ such that } V^\top C V = I}
#' where \eqn{V} is a projection matrix, i.e., \eqn{V^\top V = I}. Trace ratio optimization 
#' is pertained to various linear dimension reduction methods. It should be noted that 
#' when \eqn{C = I}, the above problem is often reformulated as a generalized eigenvalue problem 
#' since it's an easier proxy with faster computation.
#' 
#' @param A a \eqn{(p\times p)} symmetric matrix in the numerator term.
#' @param B a \eqn{(p\times p)} symmetric matrix in the denomiator term.
#' @param C a \eqn{(p\times p)} symmetric constraint matrix. If not provided, it is set as identical matrix automatically.
#' @param dim an integer for target dimension. It can be considered as the number of loadings.
#' @param method the name of algorithm to be used. Default is \code{2003Guo}.
#' @param maxiter maximum number of iterations to be performed.
#' @param eps stopping criterion for iterative algorithms. 
#' 
#' @return a named list containing
#' \describe{
#' \item{V}{a \eqn{(p\times dim)} projection matrix.}
#' \item{tr.val}{an attained maximum scalar value.}
#' }
#' 
#' @examples 
#' ## simple test
#' #  problem setting
#' p = 5
#' mydim = 2
#' A = matrix(rnorm(p^2),nrow=p); A=A%*%t(A)
#' B = matrix(runif(p^2),nrow=p); B=B%*%t(B)
#' C = diag(p)
#' 
#' #  approximate solution via determinant ratio problem formulation
#' eigAB  = eigen(solve(B,A)) 
#' V      = eigAB$vectors[,1:mydim]
#' eigval = sum(diag(t(V)%*%A%*%V))/sum(diag(t(V)%*%B%*%V))
#' 
#' #  solve using 4 algorithms
#' m12 = trio(A,B,dim=mydim, method="2012Ngo")
#' m09 = trio(A,B,dim=mydim, method="2009Jia")
#' m07 = trio(A,B,dim=mydim, method="2007Wang")
#' m03 = trio(A,B,dim=mydim, method="2003Guo")
#' 
#' #  print the results
#' line1 = '* Evaluation of the cost function'
#' line2 = paste("* approx. via determinant : ",eigval,sep="")
#' line3 = paste("* trio by 2012Ngo         : ",m12$tr.val, sep="")
#' line4 = paste("* trio by 2009Jia         : ",m09$tr.val, sep="")
#' line5 = paste("* trio by 2007Wang        : ",m07$tr.val, sep="")
#' line6 = paste("* trio by 2003Guo         : ",m03$tr.val, sep="")
#' cat(line1,"\n",line2,"\n",line3,"\n",line4,"\n",line5,"\n",line6)
#' 
#' @references 
#' \insertRef{guo_generalized_2003}{maotai}
#' 
#' \insertRef{wang_trace_2007}{maotai}
#' 
#' \insertRef{yangqingjia_trace_2009}{maotai}
#' 
#' \insertRef{ngo_trace_2012}{maotai}
#' 
#' @export
trio <- function(A, B, C, dim=2, method=c("2003Guo","2007Wang","2009Jia","2012Ngo"), maxiter=1000, eps=1e-10){
  ###################################################################
  # not completed yet.
  if (missing(C)){
    C = diag(nrow(A))
    myflag = TRUE
  } else {
    myflag = FALSE
  }
  
  if (!check_symm(A)){
    stop("* trio : an input matrix 'A' should be a square, symmetric matrix.")
  }
  if (!check_symm(B)){
    stop("* trio : an input matrix 'B' should be a square, symmetric matrix.")
  }
  if (!check_symm(C)){
    stop("* trio : an input matrix 'C' should be a square, symmetric matrix.")
  }
  sizes = rep(0,3)
  sizes[1] = nrow(A)
  sizes[2] = nrow(B)
  sizes[3] = nrow(C)
  
  if (length(unique(sizes))!=1){
    stop("* trio : all input matrices should be of same size.")
  }
  
  if (!myflag){
    eigC  = eigen(C)
    Cinv2 = eigC$vectors%*%diag(1/sqrt(eigC$values))%*%t(eigC$vectors)
    A = Cinv2%*%A%*%Cinv2
    B = Cinv2%*%B%*%Cinv2  
  }
  
  # 2009 Jia's note : B should have rank >= (m-d)
  if (as.integer(Matrix::rankMatrix(B))<(nrow(B)-dim)){
    warning("* trio : null space of 'B' is excessive. trace ratio value may diverge.")
  }
  
  
  ###################################################################
  # switch case
  V = switch(method,
             "2007Wang" = trio2007Wang(A, B, dim, eps, maxiter),
             "2003Guo"  = trio2003Guo(A, B, dim, eps, maxiter),
             "2009Jia"  = trio2009Jia(A, B, dim, eps, maxiter),
             "2012Ngo"  = trio2012Ngo(A, B, dim, eps, maxiter))
  output = list()
  output$V = V
  output$tr.val = sum(diag(t(V)%*%A%*%V))/sum(diag(t(V)%*%B%*%V))
  return(output)
}



# subroutines -------------------------------------------------------------
#' 2003 Guo et al.
#' Title : A generalized Foley-Sammon transform based on generalized fisher discriminant ...
#' 
#' @keywords internal
#' @noRd
trio2003Guo <- function(A, B, dim, eps, maxiter){
  ## translate into the language
  d = dim
  Sp = A
  Sl = B
  
  ## bisection 
  #   1. initialization
  lbd1 = 0; f1 = evalGuoDiff(lbd1, Sp, Sl, d)
  lbd2 = 1; f2 = evalGuoDiff(lbd2, Sp, Sl, d)
  if (f2 >= 0){
    while (f2 > 0){
      lbd1 = lbd2; f1 = f2;
      lbd2 = lbd2*2; f2 = evalGuoDiff(lbd2, Sp, Sl, d)  
    }
  }
  
  for (i in 1:maxiter){
    lbdm = (lbd1+lbd2)/2
    fm   = evalGuoDiff(lbdm, Sp, Sl, d)
    if (fm > 0){
      lbd1 = lbdm
      f1   = fm
    } else {
      lbd2 = lbdm
      f2   = fm
    }
    if (abs(lbd1-lbd2) < eps){
      break
    }
  }
  lbdm = (lbd1+lbd2)/2
  # W = eigen(Sp-lbdm*Sl)$vectors[,1:d] ## use RSpectra for only top 'd' ones
  W = RSpectra::eigs(Sp-lbdm*Sl,d,which="LR")$vectors
  
  ## let's try to return !
  return(W)
}


#' @keywords internal
#' @noRd
evalGuoDiff <- function(lbd, A, B, dim){
  W = RSpectra::eigs(A-lbd*B,dim,which="LR")$vectors
  # W = eigen(A-lbd*B)$vectors[,1:dim] ## use RSpectra for only top 'd' ones
  return(sum(diag(t(W)%*%(A-lbd*B)%*%W)))
}

#' 2007 Wang et al.
#' Title : Trace Ratio vs. Ratio Trace for Dimensionality Reduction
#' 
#' @keywords internal
#' @noRd
trio2007Wang <- function(A, B, dim, eps, maxiter){
  ## translate into this language
  Sp = A
  St = A+B
  m = nrow(A)
  d = dim 
  
  eigSt = base::eigen(St, symmetric = TRUE)
  mm = sum(eigSt$values > 0)
  if (mm < 1){
    stop("* (A+B) has at least one nonnegative eigenvalues.")
  }
  U  = eigSt$vectors[,1:mm]
  
  ## transform into the problem of V now.
  Spu = t(U)%*%Sp%*%U
  Stu = t(U)%*%St%*%U
  Vold = qr.Q(qr(matrix(rnorm(mm*d),ncol=d)))
  
  ## main computation
  V = cppsub_2007Wang(Vold, mm, d, Spu, Stu, maxiter, eps)
  
  ## adjust back to the original problem
  W = U%*%V
  
  ## let's try to return !
  return(W)
}

#' 2009 Jia et al 
#' Title : Trace Ratio Problem Revisited (DNM)
#' 
#' @keywords internal
#' @noRd
trio2009Jia <- function(A, B, dim, eps, maxiter){
  ## translate into the language
  d = dim
  Sp = A
  Sl = B
  
  ## Decomposed Newton Method
  lbdold = 0
  for (i in 1:maxiter){
    ## 1. compute eigendecomposition
    eigS = RSpectra::eigs(Sp-lbdold*Sl,d,which="LR")
    top.val = eigS$values   # top 'd' eigenvalues
    top.vec = eigS$vectors  # top 'd' eigenvectors
    
    ## 2. lbdnew
    lbdnew = solve2009Jia(lbdold, top.val, top.vec, Sl)
    inc    = abs(lbdnew-lbdold)
    
    ## 3. updating information
    lbdold = lbdnew
    if (inc<eps){
      break
    }
  }
  
  ## now lbdold is our solution
  W = RSpectra::eigs(Sp-lbdold*Sl,d,which="LR")$vectors
  
  ## let's try to return !
  return(W)
}


#' @keywords internal
#' @noRd
solve2009Jia <- function(lbd, eval, evec, Sl){
  d = length(eval)
  primesum = 0
  for (i in 1:d){
    tgt = as.vector(evec[,i])
    primesum = primesum - sum(as.vector(Sl%*%tgt)*tgt)
  }
  sol = (lbd*primesum - sum(eval))/primesum
  return(sol)
}

#' 2012 Ngo et al 
#' Title : The Trace Ratio Optimization Problem (Newton-Lanczos Algorithm)
#' 
#' 
#' @keywords internal
#' @noRd
trio2012Ngo <- function(A, B, dim, eps, maxiter){
  ## in the the language
  n = nrow(A)
  p = dim
  
  ## prepare the initializer
  Vold = qr.Q(qr(matrix(rnorm(n*p),ncol=p)))
  rhoold = 0
  for (i in 1:maxiter){
    Vnew   = RSpectra::eigs(A-rhoold*B,p,which="LR")$vectors
    rhonew = sum(diag(t(Vnew)%*%A%*%Vnew))/sum(diag(t(Vnew)%*%B%*%Vnew))
    
    rhoinc = abs(rhonew-rhoold)
    Vold   = Vnew
    rhoold = rhonew
    
    if (rhoinc < eps){
      break
    }
  }
  
  ## let's try to return !
  return(Vold)
}


# # 
# # simple test
# p = 100
# mydim = 10
# A = matrix(rnorm(p^2),nrow=p); A=A%*%t(A)
# B = matrix(runif(p^2),nrow=p); B=B%*%t(B)
# C = diag(p)
# 
# 
# library(geigen)
# eigAB = eigen(solve(B,A)) ## geigen(B,A) # we need largest so be careful of the order
# mylist = list()
# V = eigAB$vectors[,1:mydim]; mylist$V = V
# myval = sum(diag(t(V)%*%A%*%V))/sum(diag(t(V)%*%B%*%V)); mylist$tr.val = myval
# 
# m12 = trio(A,B,dim=mydim, method="2012Ngo")
# m9 = trio(A,B,dim=mydim, method="2009Jia")
# m7 = trio(A,B,dim=mydim, method="2007Wang")
# m3 = trio(A,B,dim=mydim, method="2003Guo")
# 
#  # try rstiefel's optimization
# library(rstiefel)
# f1 = function(w){
#   trA = sum(diag(t(w)%*%A%*%w)); trB = sum(diag(t(w)%*%B%*%w));
#   return(-trA/trB)
# }
# df1 = function(w){
#   trA = sum(diag(t(w)%*%A%*%w)); trB = sum(diag(t(w)%*%B%*%w));
#   t1 = -2*(A%*%w)*trB + trA*2*(B%*%w);
#   t2 = trB^2;
#   return(t1/t2)
# }
# f2 = function(w){
#   trA = sum(diag(t(w)%*%A%*%w)); trB = sum(diag(t(w)%*%B%*%w));
#   return(trB/trA)
# }
# df2 = function(w){
#   trA = sum(diag(t(w)%*%A%*%w)); trB = sum(diag(t(w)%*%B%*%w));
#   t1 = 2*trA*(B%*%w) - 2*trB*(A%*%w);
#   t2 = trA^2;
#   return(t1/t2)
# }
# V1 <- optStiefel(f1, df1, Vinit=rustiefel(p,mydim), maxIters = 999, tol=1e-10, verbose=TRUE)
# V2 <- optStiefel(f2, df2, Vinit=rustiefel(p,mydim), maxIters = 999, tol=1e-10, verbose=TRUE)
# print(sprintf("result of (2) inverse  : %f", 1/f2(V2)))