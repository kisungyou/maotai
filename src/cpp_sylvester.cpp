#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat cpp_sylvester(arma::mat A, arma::mat B, arma::mat C){
  arma::mat solution;
  arma::syl(solution,A,B,C);
  return(solution);
}