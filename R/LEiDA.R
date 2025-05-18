#' Leading Eigenvector Dynamics Analysis
#' 
#' Compute the leading eigenvector dynamics analysis (LEiDA) of a 
#' multivariate time series as appearing in computational neuroscience.
#' 
#' @param X A \eqn{(T,N)} matrix of multivariate time series data, where \eqn{T} is the number of time points and \eqn{N} is the number of ROIs.
#' @param TR Repetition time (in secounds)
#' @param bp Bandpass filter, a vector of length 2 with the lower and upper bounds of the bandpass filter in Hz. Default is \code{c(0.01, 0.10)}.
#' @param b_ord Butterworth order, a positive integer. Default is 2.
#' 
#' @return a list containing \describe{
#' \item{V}{A \eqn{(T,N)} matrix of the leading eigenvector time series.}
#' \item{FCD_cos}{A \eqn{(T,T)} matrix of the functional connectivity dynamics (FCD) using cosine similarity.}
#' \item{FCD_cor}{A \eqn{(T,T)} matrix of the functional connectivity dynamics (FCD) using Pearson correlation.}
#' }
#' 
#' @export
LEiDA <- function(X, TR, bp=c(0.01, 0.10), b_ord = 2){
  # 1. PREPROCESSING
  fs <- 1/as.double(TR)
  par_T <- base::nrow(X)
  par_N <- base::ncol(X)
  
  Wn <- bp/(fs/2) # normalised cut-offs
  bpFilt <- gsignal::butter(round(b_ord), Wn, type="pass")
  
  # 2. DEMEANING, DETRENDING, AND FILTERING
  Xclean <- array(0,c(par_T, par_N)) # filtered data
  for (j in seq_len(par_N)){         # run for each ROI
    ts <- as.vector(X[,j])
    ts <- ts - base::mean(ts)                  # demeaning
    ts <- pracma::detrend(ts, tt = "linear")   # detrending
    Xclean[,j] = gsignal::filtfilt(bpFilt, ts) # filtering
  }
  
  # 3. PHASE 
  Phase <- array(0,c(par_T, par_N))
  for (j in seq_len(par_N)){
    Phase[,j] = base::Arg(gsignal::hilbert(Xclean[,j])) # phase
  }
  
  # 4. COMPUTE AND RETURN
  computed = src_leida(Phase)
  return(computed)
}
