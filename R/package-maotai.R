#' Tools for Matrix Algebra, Optimization and Inference
#' 
#' Matrix is an universal and sometimes primary object/unit in applied mathematics and statistics. We provide a number of algorithms for selected problems in optimization and statistical inference.
#' This package contains following functions, 
#' \tabular{ll}{
#' FUNCTION \tab DESCRIPTION \cr
#' \code{\link{lgpa}} \tab Large-scale Generalized Procrustes Analysis \cr
#' \code{\link{lyapunov}} \tab Solve Lyapunov Equation \cr
#' \code{\link{matderiv}} \tab Numerical Approximation to Gradient of a Function with Matrix Argument \cr
#' \code{\link{pdeterminant}} \tab Calculate the Pseudo-Determinant of a Matrix \cr
#' \code{\link{shortestpath}} \tab Find Shortest Path using Floyd-Warshall Algorithm \cr
#' \code{\link{sylvester}} \tab Solve Sylvester Equation \cr
#' \code{\link{trio}} \tab Trace Ratio Optimation
#' }
#'
#' @docType package
#' @name maotai
#' @aliases maotai-package
#' @import Rdpack
#' @importFrom stats rnorm
#' @importFrom shapes procGPA
#' @importFrom utils packageVersion
#' @importFrom RSpectra eigs
#' @importFrom Matrix rankMatrix
#' @importFrom Rcpp evalCpp
#' @useDynLib maotai
NULL