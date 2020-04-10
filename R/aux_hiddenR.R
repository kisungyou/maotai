# Hidden Functions for Future Use
# these functions can be loaded using 'utils::getFromNamespace'
# by the command 'getFromNamespace("function_name","maotai");
# 
# 01. hidden_pinv        : pseudo-inverse
# 02. hidden_vech        : half vectorization including the diagonal.
#     hidden_ivech         diagonal inclusion is also taken care. 
# 03. hidden_lab2ind     : create an index list from a label vector
#     hidden_ind2lab       given an index list, create a label vector


# 01. hidden_pinv ---------------------------------------------------------
#' @keywords internal
hidden_pinv <- function(A){
  return(aux_pinv(A))
}

# 02. hidden_vech & hidden_ivech ------------------------------------------
#' @keywords internal
hidden_vech <- function(A, diag=TRUE){
  if ((!is.matrix(A))||(nrow(A)!=ncol(A))){
    stop("* hidden_vech : input should be a square matrix.")
  }
  mydiag = as.logical(diag)
  return(A[base::lower.tri(A, diag=mydiag)])
}
#' @keywords internal
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
hidden_ind2lab <- function(index){
  K = length(index)
  N = sum(unlist(lapply(index, length)))
  
  output = rep(0,N)
  for (k in 1:K){
    output[index[[k]]] = k
  }
  return(output)
}