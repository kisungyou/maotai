# Hidden Functions for Future Use
# these functions can be loaded using 'utils::getFromNamespace'
# by the command 'getFromNamespace("function_name","maotai");
# 
# 01. hidden_pinv         : pseudo-inverse
# 02. hidden_vech         : half vectorization including the diagonal.
#     hidden_ivech          diagonal inclusion is also taken care. 
# 03. hidden_lab2ind      : create an index list from a label vector
#     hidden_ind2lab        given an index list, create a label vector
# 04. hidden_subsetid     : generate split of the subset id
# 05. hidden_geigen       : do 'geigen' operation; Decreasing order
# 06. hidden_knn
# 07. hidden_knee_clamped : knee-point detection with clamped least squares - return idx
# 08. hidden_knn_binary   : return a sparse binary matrix for Euclidean data excluding selfl dgCMatrix


# 01. hidden_pinv ---------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_pinv <- function(A){
  return(aux_pinv(A))
}

# 02. hidden_vech & hidden_ivech ------------------------------------------
#' @keywords internal
#' @noRd
hidden_vech <- function(A, diag=TRUE){
  if ((!is.matrix(A))||(nrow(A)!=ncol(A))){
    stop("* hidden_vech : input should be a square matrix.")
  }
  mydiag = as.logical(diag)
  return(A[base::lower.tri(A, diag=mydiag)])
}
#' @keywords internal
#' @noRd
hidden_ivech <- function(a, diag=TRUE){
  k = length(a)
  if (diag){
    n = round((-1 + sqrt(1+(8*k)))/2)
    output = array(0,c(n,n))
    output[lower.tri(output, diag = TRUE)] = a
    output = output + t(output)
    diag(output) = diag(output)/2
  } else {
    n = round((1+sqrt(1+8*k))/2)
    output = array(0,c(n,n))
    output[lower.tri(output, diag = FALSE)] = a
    output = output + t(output)
  }
  return(output)
}


# 03. hidden_lab2ind & hidden_ind2lab -------------------------------------
#' @keywords internal
#' @noRd
hidden_lab2ind <- function(label){
  ulabel = base::unique(label)
  nlabel = length(ulabel)
  
  index  = list()
  for (k in 1:nlabel){
    index[[k]] = which(label==ulabel[k])
  }
  return(index)
}
#' @keywords internal
#' @noRd
hidden_ind2lab <- function(index){
  K = length(index)
  N = sum(unlist(lapply(index, length)))
  
  output = rep(0,N)
  for (k in 1:K){
    output[index[[k]]] = k
  }
  return(output)
}

# 04. hidden_subsetid -----------------------------------------------------
#' @keywords internal
#' @noRd
hidden_subsetid <- function(n, k){
  return(base::split(base::sample(n), base::sort(n%%k)))
}

# 05. hidden_geigen -------------------------------------------------------
#' It mimics the behavior of 'geigen' function with normalization added
#' @keywords internal
#' @noRd
hidden_geigen <- function(A, B, normalize=TRUE){
  n    = nrow(A)
  runs = cpp_geigen(A,B)
  
  tval = as.vector(base::Re(runs$values))[n:1]
  tvec = base::Re(runs$vectors)[,n:1]
  if (normalize){
    for (i in 1:n){
      tgt = as.vector(tvec[,i])
      tvec[,i] = tgt/sqrt(sum(tgt^2))
    }  
  }

  output = list()
  output$values  = tval
  output$vectors = tvec
  return(output)
}

# 06. hidden_knn ----------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_knn <- function(dat, nnbd=2, ...){
  nnbd = round(nnbd)
  return(RANN::nn2(dat, k=nnbd, ...))
}


# 07. hidden_knee_clamped -------------------------------------------------
#' @keywords internal
#' @noRd
hidden_knee_clamped_basic <- function(x, y){
  m = length(x)
  c = x[1]
  d = y[1]
  a = x[m]
  b = y[m]
  
  y2 = (((b-d)/(a-c))*(x-c))+d
  return(sum((y-y2)^2))
}
#' @keywords internal
#' @noRd
hidden_knee_clamped <- function(x, y){
  x = as.vector(x)
  y = as.vector(y)
  n = length(x)
  if (n < 3){
    stop("* knee point detection : length must be larger than 2.")
  }
  scores = rep(Inf, n)
  for (i in 2:(n-1)){
    x.left = x[1:i]
    y.left = y[1:i]
    
    x.right = x[i:n]
    y.right = y[i:n]
    
    term1 = hidden_knee_clamped_basic(x.left, y.left)
    term2 = hidden_knee_clamped_basic(x.right, y.right)
    scores[i] = term1+term2
  }
  return(which.min(scores)) # return the index of the minimal SSE's
}

# 08. hidden_knn_binary ---------------------------------------------------
#     excluding self and return a binary sparse adjacency matrix
#' @keywords internal
#' @noRd
hidden_knn_binary <- function(data, nnbd=1, construction=c("or","and")){
  n = base::nrow(data)
  nnbd = max(round(nnbd), 1)
  construction = match.arg(construction)
  if (all(construction=="and")){
    intersect = TRUE
  } else {
    intersect = FALSE
  }
  
  run_knn = RANN::nn2(data, k=nnbd+1)$nn.idx[,2:(nnbd+1)]-1 # -1 for C++ convention
  run_res = src_construct_by_knn(run_knn, intersect)
  return(run_res)
}
