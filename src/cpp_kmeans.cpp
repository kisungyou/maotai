#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List cpp_kmeans(arma::mat data, int k){ // rows are stacked observations
  // parameters
  int N = data.n_rows;
  int p = data.n_cols;
  int csub = 10; // cardinality of a random subset
  if (N/2 < 10){
    csub = N/2;
  }
  
  // prepare for kmeans
  arma::mat means(p,k,fill::zeros); // armadillo reference
  bool status = arma::kmeans(means, data.t(), k, random_subset, csub, false);
  if (status==false){
    Rcpp::stop("* epmeans : k-means failed.");
  }
  return Rcpp::List::create(Rcpp::Named("means")=means.t());
}