#' Install 'SciPy' Python Module
#' 
#' 
#' @keywords internal
#' @noRd
install_scipy <- function(method = "auto", conda = "auto") {
  reticulate::py_install("scipy", method = method, conda = conda)
}