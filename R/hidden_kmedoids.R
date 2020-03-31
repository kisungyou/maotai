#' Hidden Function for 'k-medoids'
#' 
#' 
#' @keywords internal
#' @noRd
hidden_kmedoids <- function(distobj, nclust=2){
  myk = round(nclust)
  return(cluster::pam(distobj, k = myk))
}