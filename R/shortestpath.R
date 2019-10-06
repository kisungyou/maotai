#' Find Shortest Path using Floyd-Warshall Algorithm
#'
#' This is a fast implementation of Floyd-Warshall algorithm to find the
#' shortest path in a pairwise sense using \code{RcppArmadillo}. A logical input
#' is also accepted. The given matrix should contain pairwise distance values \eqn{d_{i,j}} where 
#' \eqn{0} means there exists no path for node \eqn{i} and {j}.
#'
#' @param dist either an \eqn{(n\times n)} matrix or a \code{dist} class object.
#' 
#' @return an \eqn{(n\times n)} matrix containing pairwise shortest path length.
#' 
#' @examples
#' ## simple example : a ring graph
#' #  edges exist for pairs
#' A = array(0,c(10,10))
#' for (i in 1:9){
#'   A[i,i+1] = 1
#'   A[i+1,i] = 1
#' }
#' A[10,1] <- A[1,10] <- 1
#' 
#' # compute shortest-path and show the matrix
#' sdA <- shortestpath(A)
#' image(sdA, main="shortest path length for ring graph")
#'
#' @references 
#' \insertRef{floyd_algorithm_1962}{maotai}
#' 
#' \insertRef{warshall_theorem_1962}{maotai}
#' 
#' @export
shortestpath <- function(dist){
  input = dist
  # class determination
  if (class(dist)=="matrix"){
    distnaive = dist
  } else if (class(dist)=="dist"){
    distnaive = as.matrix(dist)
  } else {
    stop("* shortestpath : input 'dist' should be either (n*n) matrix or 'dist' class object.")
  }
  # consider logical input
  if (any(is.logical(distnaive))){
    distnaive = distnaive*1
  }
  # set as -Inf for 0 values
  mepsil  = .Machine$double.eps
  distnaive[which(distnaive<5*mepsil)] = -Inf
  distgeo   = aux_shortestpath(distnaive)
  return(distgeo)
}