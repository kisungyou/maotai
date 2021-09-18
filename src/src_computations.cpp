#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/* Functions for Generic Computations
 * 
 * (01) src_construct_by_knn [hidden_knn_binary]
 * 
 * 
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