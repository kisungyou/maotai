#' Find a median from the pairwise dissimilarity
#' 
#' 
#' @keywords internal
#' @noRd
findmed <- function(d, method=c("geometric","metric")){
  ######################################################
  # Preprocessing
  if (inherits(d, "dist")){
    d = as.matrix(d)
  } else {
    if (!is.matrix(d)){
      stop("* findmed : input 'd' should be a matrix.")
    }
  }
  if (missing(method)){
    mymethod = "geometric"
  } else {
    mymethod = match.arg(method)
  }
  diag(d) = 0
  
  ######################################################
  # Compute and Return
  if (all(mymethod=="metric")){ # metric median;
    return(which.min(base::rowSums(d)))
  } else {                      # geometric median;
    m = base::nrow(d)
    vecm = rep(0,m)
    for (i in 1:m){
      tgt = base::sort(as.vector(d[i,]), decreasing=FALSE)
    }
  }
}