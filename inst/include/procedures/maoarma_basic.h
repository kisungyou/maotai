#ifndef MAOTAI_MAOARMA_BASIC_H
#define MAOTAI_MAOARMA_BASIC_H

#define ARMA_NO_DEBUG

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include "../include/maoarma.hpp"
#include <cassert>

using namespace Rcpp;
using namespace arma;
using namespace std;

// auxiliary functiona
// * aux_pnorm : given a row vector, compute p-norm
inline double aux_pnorm(arma::rowvec x, double p){
  assert(p>0);
  int N = x.n_elem;
  double output = 0.0;
  for (int n=0;n<N;n++){
    output += std::pow(static_cast<float>(x(n)), p);
  }
  return(std::pow(static_cast<float>(output), 1.0/p));
}
//--------------------------------------------------------//
// main namespace functions
inline arma::mat maospace::pdist(arma::mat X, double p){
  // prepare
  int N = X.n_rows;
  arma::mat output(N,N,fill::zeros);
  // iterate
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      output(i,j) = aux_pnorm(X.row(i)-X.row(j),p);
      output(j,i) = output(i,j);
    }
  }
  // return
  return(output);
}
inline arma::mat maospace::pdist2(arma::mat X, arma::mat Y, double p){
  // prepare
  int N = X.n_rows;
  int M = Y.n_rows;
  arma::mat output(N,M,fill::zeros);
  // iterate
  for (int n=0;n<N;n++){
    for (int m=0;m<M;m++){
      output(n,m) = aux_pnorm(X.row(n)-Y.row(m), p);
    }
  }
  // return
  return(output);
}
inline arma::umat maospace::knn(arma::mat X, int k, double p){
  // parameter
  assert(k>0);
  int n = X.n_rows;
  arma::umat indices(n,k,fill::zeros);
  arma::vec distvec(n,fill::zeros);
  arma::uvec sorted;
  
  for (int i=0;i<n;i++){
    distvec.fill(0.0);
    for (int j=0;j<n;j++){
      distvec(j) = aux_pnorm(X.row(i)-X.row(j),p);
    }
    sorted = arma::sort_index(distvec,"ascend");
    indices.row(i) = arma::trans(sorted.rows(1,k));
  }
  return(indices);
}
inline double maospace::trace(arma::mat X){
  // parameter
  int n = X.n_rows;
  int m = X.n_cols;
  assert(n==m);
  
  double output = 0.0;
  for (int i=0;i<n;i++){
    output += X(i,i);
  }
  return(output);
}

#endif