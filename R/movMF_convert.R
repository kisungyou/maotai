#' Convert 'movMF' object
#' 
#' Given an output from the movMF package's movMF function, 
#' convert them into the standard mixture parameter format.
#' 
#' @param movMF_object a movMF object of \eqn{K} components in \eqn{d} dimensions.
#'
#' @return a named list containing \describe{
#' \item{means}{a \eqn{(K \times d)} matrix of means}
#' \item{concentrations}{a \eqn{K} vector of concentrations}
#' \item{weights}{a \eqn{K} vector of weights}
#' }
#' 
#' @export
movMF_convert <- function(movMF_object){
  ###############################################
  # Preprocessing
  if (!inherits(movMF_object, "movMF")){
    stop("* movMF_convert : Input object is not a 'movMF' object")
  }
  
  ###############################################
  # Change
  output = list()
  output$means <- movMF_object$theta/sqrt(rowSums(movMF_object$theta^2))
  output$concentrations <- sqrt(rowSums(movMF_object$theta^2))
  output$weights <- movMF_object$alpha
  return(output)
}
