#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat solve_lyapunov(arma::mat A, arma::mat B, arma::mat C){
  // simply solve it !
  arma::mat solution;
  arma::syl(solution, A, B, C);
  return(solution);
}
