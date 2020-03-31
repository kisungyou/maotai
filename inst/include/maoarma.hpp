// #ifndef MAOTAI_MAOARMA_HPP
// #define MAOTAI_MAOARMA_HPP
// 
// #pragma once

#include <RcppArmadillo.h>
#include "procedures/maoarma_basic.h"

// available functions
//  (basic:umat)   knn       : return (n x k) index vector. 0 starting C convention.
//  (basic:mat)    pdist     : return (n x n) distance matrix
//  (basic:mat)    pdist2    : return (n x m) distance matrix for X(n) and Y(m)
//  (basic:double) trace     : given a square matrix, compute the trace

// auxiliary functions
//  (basic:double) aux_pnorm : compute p-norm

namespace maospace{
  arma::umat knn(arma::mat X, int k, double p);
  arma::mat pdist(arma::mat X, double p);
  arma::mat pdist2(arma::mat X, arma::mat Y, double p);
  double trace(arma::mat X);
}

// How to use this 'hpp' file 
// * set 'import' and 'linkingto' with maotai
// * In C++ scripts, do the followings
//    - #include "maoarma.hpp"
//    - // [[Rcpp::depends(maotai)]]
//    - using namespace maospace

// #endif