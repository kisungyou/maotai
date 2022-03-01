#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * eval_gaussian      : evaluate single   observation
 * eval_gaussian_data : evaluate multiple observations 
 * eval_gmm           : evaluate mixture of a single
 * eval_gmm_data      : evaluate mixture of multiple
 */


// [[Rcpp::export]]
double eval_gaussian(arma::vec x, arma::vec mu, arma::mat cov){
  // some constants
  double pi2 = 6.28318530717958623199592693708837032318115234375;
  double d = static_cast<double>(x.n_elem);
  
  arma::vec xdiff = x-mu;
  arma::vec xsolv = arma::solve(cov, xdiff);
  
  double term1 = arma::dot(xdiff,xsolv)*(-0.5);
  double term2 = -(d/2.0)*std::log(pi2);
  double term3 = -0.5*std::log(static_cast<float>(arma::det(cov)));
  
  return(std::exp(term1+term2+term3));
}

// [[Rcpp::export]]
arma::vec eval_gaussian_data(arma::mat X, arma::vec mu, arma::mat cov){ // (n x p) convention
  // some constants
  double pi2 = 6.28318530717958623199592693708837032318115234375;
  int n = X.n_rows;
  int d = mu.n_elem;

  // output
  arma::vec output(n,fill::zeros);
  
  // common objects
  double term1 = 0.0;
  double term2 = -(static_cast<double>(d)/2.0)*std::log(pi2);
  double term3 = -0.5*std::log(static_cast<float>(arma::det(cov)));

  arma::rowvec xtgtrow(d,fill::zeros);
  arma::colvec xtgtcol(d,fill::zeros);
  arma::vec xdiff(d,fill::zeros);
  arma::vec xsolv(d,fill::zeros);
  
  for (int i=0;i<n;i++){
    xtgtrow = X.row(i);
    xtgtcol = X.row(i).t();
    xdiff = xtgtcol-mu;
    xsolv = arma::solve(cov, xdiff);
    term1 = arma::dot(xdiff,xsolv)*(-0.5);
    output(i) = std::exp(term1 + term2 + term3);
  }
  return(output);
}

/* X      : (N x p) data observation
 * mus    : (K x p) means
 * covs   : (p x p x K) covariances
 * weight : (K) vector of weights
 */
// [[Rcpp::export]]
arma::vec eval_gmm_data(arma::mat X, arma::mat mus, arma::cube covs, arma::vec weight){
  // parameters
  int N = X.n_rows;
  int p = X.n_cols;
  int K = weight.n_elem;
  
  // iterate over different clusters
  arma::vec mu(p,fill::zeros);
  arma::mat cov(p,p,fill::zeros);
  arma::vec tmpvals(N,fill::zeros);
  arma::mat records(N,K,fill::zeros);
  for (int k=0;k<K;k++){
    mu  = mus.row(k).t(); // per-class mean
    cov = covs.slice(k);  // per-class covariance
    tmpvals = eval_gaussian_data(X, mu, cov); // from a single density
    records.col(k) = tmpvals*weight(k);
  }
  
  // return a rowsum
  arma::vec output = arma::sum(records, 1);
  return(output);
}

// [[Rcpp::export]]
double eval_gmm(arma::vec x, arma::mat mus, arma::cube covs, arma::vec weight){
  int p = x.n_elem;
  arma::mat X = arma::reshape(x, 1, p);
  return(arma::accu(eval_gmm_data(X, mus, covs, weight)));
}
