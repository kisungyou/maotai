# CHECKERS ----------------------------------------------------------------
# 01. check_sqmat : check if a square matrix
# 02. check_symm  : check if a square, symmetric matrix




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