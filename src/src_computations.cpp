#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* Functions for Generic Computations
 * 
 * (01) src_construct_by_knn [hidden_knn_binary]
 * (02) src_gaussbary_2002R : Barycenter of Gaussian Covariances
 * (03) src_gaussbary_2016A : Barycenter of Gaussian Covariances
 * (04) src_cov2corr        : use Wikipedia's formula
 */

// (01) src_construct_by_knn [hidden_knn_binary] ===============================
// [[Rcpp::export]]
arma::sp_umat src_construct_by_knn(arma::umat& nn_idx, bool intersection){
  // prepare
  int n = nn_idx.n_rows;
  int k = nn_idx.n_cols;
  arma::sp_umat output(n,n);
  
  // fill in 
  for (int i=0; i<n; i++){
    for (int j=0; j<k; j++){
      output(i,nn_idx(i,j)) = 1;
    }
  }
  
  if (intersection){ // case : intersection 
    for (int i=0; i<(n-1); i++){
      for (int j=(i+1); j<n; j++){
        if ((output(i,j) > 0)&&(output(j,i) > 0)){
          output(i,j) = 1;
          output(j,i) = 1;
        } else {
          output(i,j) = 0;
          output(j,i) = 0;
        }
      }
    }
  } else {           // case : union
    for (int i=0; i<(n-1); i++){
      for (int j=(i+1); j<n; j++){
        if ((output(i,j) > 0)||(output(j,i) > 0)){
          output(i,j) = 1;
          output(j,i) = 1;
        } else {
          output(i,j) = 0;
          output(j,i) = 0;
        }
      }
    }
  }
  return(output);
}

// (02) src_gaussbary_2002R ====================================================
// [[Rcpp::export]]
Rcpp::List src_gaussbary_2002R(arma::cube &array3d, arma::vec &weight, int maxiter, double abstol){
  // PREPARE
  int p = array3d.n_rows;
  int N = array3d.n_slices;
  
  double S_inc = 10000.0;
  arma::mat S_old = arma::mean(array3d, 2);
  int S_old_rank = arma::rank(S_old);
  if (S_old_rank < p){
    S_old.fill(0.0);
    for (int n=0; n<N; n++){
      S_old = arma::logmat_sympd(array3d.slice(n))/static_cast<double>(N);
    }
    S_old = arma::expmat_sym(S_old);
  }
  
  arma::mat S_new(p,p,fill::zeros);
  arma::mat S_oldhalf(p,p,fill::zeros);
  
  arma::vec S_weight = weight/arma::accu(weight);
  
  // MAIN
  int itcount = 0;
  for (int it=0; it<maxiter; it++){
    // MAIN-1. initialize 'new' and apply 'sqrtm' on the old one.
    itcount += 1;
    S_new.fill(0.0);
    S_oldhalf = arma::sqrtmat_sympd(S_old);
    
    // MAIN-2. iteratively update
    for (int n=0; n<N; n++){
      S_new += S_weight(n)*arma::sqrtmat_sympd(S_oldhalf*array3d.slice(n)*S_oldhalf);
    }
    
    // MAIN-3. update
    S_inc = arma::norm(S_old-S_new,"fro");
    S_old = S_new;
    if (S_inc < abstol){
      break;
    }
  }
  
  // RETURN
  return(Rcpp::List::create(Rcpp::Named("mean")=S_old,
                            Rcpp::Named("iter")=itcount));
}

// (03) src_gaussbary_2016A ====================================================
// [[Rcpp::export]]
Rcpp::List src_gaussbary_2016A(arma::cube &array3d, arma::vec &weight, int maxiter, double abstol){
  // PREPARE
  int p = array3d.n_rows;
  int N = array3d.n_slices;
  
  double S_inc = 10000.0;
  arma::mat S_old = arma::mean(array3d, 2);
  int S_old_rank = arma::rank(S_old);
  if (S_old_rank < p){
    S_old.fill(0.0);
    for (int n=0; n<N; n++){
      S_old = arma::logmat_sympd(array3d.slice(n))/static_cast<double>(N);
    }
    S_old = arma::expmat_sym(S_old);
  }

  arma::mat S_tmp(p,p,fill::zeros);
  arma::mat S_new(p,p,fill::zeros);
  
  arma::mat S0_half(p,p,fill::zeros);
  arma::mat S0_hinv(p,p,fill::zeros);
  
  arma::vec S_weight = weight/arma::accu(weight);
  
  // MAIN
  int itcount = 0;
  for (int it=0; it<maxiter; it++){
    // MAIN-1. initialize 'tmp' and compute auxiliary ones
    itcount += 1;
    S_tmp.fill(0.0);
    S0_half = arma::sqrtmat_sympd(S_old);
    S0_hinv = arma::inv_sympd(S0_half);
    
    // MAIN-2. compute the intermediate summation
    for (int n=0; n<N; n++){
      S_tmp += S_weight(n)*arma::sqrtmat_sympd(S0_half*array3d.slice(n)*S0_half);
    }
    
    // MAIN-3. finally attain an iterate
    S_new = S0_hinv*S_tmp*S_tmp*S0_hinv;
    
    // MAIN-3. update
    S_inc = arma::norm(S_old-S_new,"fro");
    S_old = S_new;
    if (S_inc < abstol){
      break;
    }
  }
  
  // RETURN
  return(Rcpp::List::create(Rcpp::Named("mean")=S_old,
                            Rcpp::Named("iter")=itcount));
}

// (04) src_cov2corr ===========================================================
// [[Rcpp::export]]
arma::mat src_cov2corr(arma::mat &covmat){
  int N = covmat.n_rows;
  arma::mat precision = arma::inv_sympd(covmat);
  arma::mat output(N,N,fill::ones);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output(i,j) = precision(i,j)/std::sqrt(precision(i,i)*precision(j,j));
      output(j,i) = output(i,j);
    }
  }
  return(output);
}
