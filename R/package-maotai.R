#' Tools for Matrix Algebra, Optimization and Inference
#' 
#' A matrix is an universal and sometimes primary object/unit in applied mathematics and statistics. We provide a number of algorithms for selected problems in optimization and statistical inference.
#' 
#' @keywords internal
#' @name package-maotai
#' @import Rdpack
#' @noRd
#' @importFrom dbscan dbscan
#' @importFrom fastcluster hclust
#' @importFrom RANN nn2
#' @importFrom cluster pam silhouette agnes
#' @importFrom stats as.dist knots ecdf rnorm runif quantile dist rgamma rgeom var cov lm
#' @importFrom shapes procGPA
#' @importFrom Rtsne Rtsne
#' @importFrom pracma cross detrend
#' @importFrom utils packageVersion
#' @importFrom RSpectra eigs
#' @importFrom Matrix rankMatrix
#' @importFrom Rcpp evalCpp
#' @importFrom gsignal hilbert butter filtfilt
#' @useDynLib maotai
"_PACKAGE"
# pack <- "maotai"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))
