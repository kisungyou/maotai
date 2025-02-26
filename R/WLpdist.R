#' Pairwise Wasserstein-like Distance between two vMF distributions
#' 
#' Given a collection of von Misees-Fisher (vMF) distributions, compute the pairwise 
#' distance using the Wasserstein-like distance from an approximate Wasserstein geometry.
#' 
#' @param means An \eqn{(n \times p)} matrix where each row represents the mean 
#' direction of one of the \eqn{n} vMF distributions.
#' @param concentrations A length-\eqn{n} vector of nonnegative concentration parameters.
#' 
#' @return An \eqn{(n \times n)} matrix of pairwise distances.
#' 
#' @examples
#' \donttest{
#' # Set seed for reproducibility
#' set.seed(123)
#' 
#' # Generate two classes of mean directions around north and south poles
#' means1 = array(0,c(50,2)); means1[,2] = rnorm(50, mean=1, sd=0.25)
#' means2 = array(0,c(50,2)); means2[,2] = rnorm(50, mean=-1, sd=0.25)
#' means1 = means1/sqrt(rowSums(means1^2))
#' means2 = means2/sqrt(rowSums(means2^2))
#' 
#' # Concatenate the mean directions
#' data_means = rbind(means1, means2)
#' 
#' # Generate concentration parameters
#' data_concentrations = rnorm(100, mean=20, sd=1)
#' 
#' # Compute the pairwise distance matrix
#' pdmat = WLpdist(data_means, data_concentrations)
#' 
#' # Visualise the pairwise distance matrix
#' opar <- par(no.readonly=TRUE)
#' image(pdmat, main="Pairwise Wasserstein-like Distance")
#' par(opar)
#' }
#' 
#' @export
WLpdist <- function(means, concentrations){
  ###############################################
  # Preprocessing
  # means
  if (!is.matrix(means)){
    data_means = cpp_WL_normalise(as.matrix(means))
  } else {
    data_means = cpp_WL_normalise(means)
  }
  
  # concentrations
  data_kappa = as.vector(concentrations)
  if (length(data_kappa)!=base::nrow(data_means)){
    stop("* WLpdist : cardinalities of the means and concentrations do not match.")
  }
  
  # parameters
  N = base::nrow(means)
  d = base::ncol(means)
  
  
  ###############################################
  # COMPUTATION
  output = array(0,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      term1 = sphere_dist(as.vector(data_means[i,]), as.vector(data_means[j,]))
      term2sq = (d-1)*((1/sqrt(data_kappa[i]) - 1/sqrt(data_kappa[j]))^2)
      output[i,j] <- output[j,i] <- sqrt((term1*term1) + term2sq)
    }
  }
  
  ###############################################
  # RETURN
  return(output)
}
