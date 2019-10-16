#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec dat2centers(arma::rowvec data, arma::mat &centers){
  // parameters
  int K = centers.n_rows;
  int p = data.n_cols;
  
  // compute
  arma::vec dic(K,fill::zeros);
  arma::rowvec diffvec(p,fill::zeros);
  for (int k=0;k<K;k++){
    diffvec = data-centers.row(k);
    dic(k)  = arma::dot(diffvec,diffvec);
  }
  
  // report
  return(dic);
}
