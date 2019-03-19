# CHECKERS ----------------------------------------------------------------
# 01. check_sqmat : check if a square matrix




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

