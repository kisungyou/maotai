#' Cayley-Menger Determinant
#' 
#' Cayley-Menger determinant is a formula of a \eqn{n}-dimensional simplex 
#' with respect to the squares of all pairwise distances of its vertices.
#' 
#' @param data an \eqn{(n\times p)} matrix of row-stacked observations.
#' @return a list containing\describe{
#' \item{det}{determinant value.}
#' \item{vol}{volume attained from the determinant.}
#' }
#' 
#' @examples 
#' ## USE 'IRIS' DATASET
#' data(iris)
#' X = as.matrix(iris[,1:4])
#' 
#' ## COMPUTE CAYLEY-MENGER DETERMINANT
#' #  since k=4 < n=149, it should be zero.
#' cayleymenger(X)
#' 
#' @export
cayleymenger <- function(data){
  # Preprocessing 
  if (!check_datamat(data)){
    stop("* cayleymenger : an input 'data' should be a matrix without any missing/infinite values.")
  }

  # compute pairwise distance
  Dtmp = stats::as.dist(cpp_pdist(data))
  
  # compute and return
  return(cayleymenger_internal(Dtmp))
}

#' @keywords internal
#' @noRd
cayleymenger_internal <- function(distobj){
  Dold = (as.matrix(distobj)^2)
  n    = base::nrow(Dold)
  Dnew = rbind(cbind(Dold, rep(1,n)), c(rep(1,n),0))
  
  val.det = base::det(Dnew)
  n       = n+1
  val.vol = base::sqrt(base::exp(base::log(((-1)^(n+1))*val.det) - n*log(2) - (2*base::lfactorial(n))))
  
  output = list(det=val.det, vol=val.vol)
  return(output)
}