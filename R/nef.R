#' Negative Eigenfraction
#' 
#' Negative Eigenfraction (NEF) is a measure of distortion for the data 
#' whether they are lying in Euclidean manner or not. When the value is exactly 0, it means 
#' the data is Euclidean. On the other hand, when NEF is far away from 0, it means not Euclidean. 
#' The concept of NEF is closely related to the definiteness of a Gram matrix. 
#' 
#' @param data an \eqn{(n\times p)} matrix whose rows are observations.
#' 
#' @return a nonnegative NEF value.
#' 
#' @examples 
#' ## use simple example of iris dataset 
#' data(iris) 
#' mydat = as.matrix(iris[,1:4])
#' 
#' ## calculate NEF
#' nef(mydat)
#' 
#' @references 
#' \insertRef{pekalska_noneuclidean_2006}{maotai}
#' 
#' @export
nef <- function(data){
  ############################################################
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* nef : an input 'data' should be a matrix without any missing/infinite values.")
  }
  xdiss = stats::as.dist(cpp_pdist(data))
  
  ############################################################
  # Compute and Return
  output = hidden_nef(xdiss)
  return(output)
}