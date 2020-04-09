#' Negative Eigenvalue Magnitude
#' 
#' Negative Eigenvalue Magnitude (NEM) is a measure of distortion for the data 
#' whether they are lying in Euclidean manner or not. When the value is exactly 0, it means 
#' the data is Euclidean. On the other hand, when NEM is far away from 0, it means not Euclidean. 
#' The concept of NEM is closely related to the definiteness of a Gram matrix. 
#' 
#' @param data an \eqn{(n\times p)} matrix whose rows are observations.
#' 
#' @return a nonnegative NEM value.
#' 
#' @examples 
#' ## use simple example of iris dataset 
#' data(iris) 
#' mydat = as.matrix(iris[,1:4])
#' 
#' ## calculate NEM
#' nem(mydat)
#' 
#' @references 
#' \insertRef{pekalska_noneuclidean_2006}{maotai}
#' 
#' @export
nem <- function(data){
  ############################################################
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* nem : an input 'data' should be a matrix without any missing/infinite values.")
  }
  xdiss = stats::as.dist(cpp_pdist(data))
  
  ############################################################
  # Compute and Return
  output = hidden_nem(xdiss)
  return(output)
}
