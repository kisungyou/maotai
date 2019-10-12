# CHECKERS ----------------------------------------------------------------
# 01. check_sqmat    : if a square matrix
# 02. check_symm     : if a square, symmetric matrix
# 03. check_datalist : if a list of same-dimensional data



# 01. check_sqmat ---------------------------------------------------------
#' @keywords internal
#' @noRd
check_sqmat <- function(x){
  cond1 = is.matrix(x)
  cond2 = (nrow(x)==ncol(x))
  cond3 = (!(any(is.infinite(x))||any(is.null(x))))
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# 02. check_symm ----------------------------------------------------------
#' @keywords internal
#' @noRd
check_symm <- function(x){
  cond1 = check_sqmat(x)
  cond2 = isSymmetric(x)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# 03. check_datalist : if a list of multivariate data

# 03. check_datalist ------------------------------------------------------
#' @keywords internal
#' @noRd
check_datalist <- function(dlist){
  cond1 = (is.list(dlist))
  if (is.vector(dlist[[1]])){
    cond2 = all(unlist(lapply(dlist, is.vector))==TRUE)
    if (cond1&&cond2){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    cond2 = all(unlist(lapply(dlist, is.matrix))==TRUE)
    cond3 = (length(unique(unlist(lapply(dlist, ncol))))==1)
    if (cond1&&cond2&&cond3){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}
