#include "RcppArmadillo.h"
#include "evaluations.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


/*
 * 1. cpp_pairwise_L2 : L2 distance between GMM's.
 * 2. integrate_1d    : 1d integration, generic but used in distance computation
 */

///////////////////////////////////////////////////////////////////
// 1. cpp_pairwise_L2
// [[Rcpp::export]]
Rcpp::List cpp_pairwise_L2(arma::mat muA, arma::mat muB, arma::cube covA, arma::cube covB){
  // parameters
  int N = muA.n_rows;
  int M = muB.n_rows;
  int p = muA.n_cols;
  
  // output
  arma::mat matA(N,N,fill::zeros);
  arma::mat matB(M,M,fill::zeros);
  arma::mat matC(N,M,fill::zeros);
  
  // preparation
  arma::vec parvec1(p,fill::zeros);
  arma::vec parvec2(p,fill::zeros);
  arma::mat parcovv(p,p,fill::zeros);
  
  // matA : N normals from mixture 1
  for (int n=0;n<N;n++){
    parvec1 = muA.row(n).t();
    parcovv = 2*covA.slice(n);
    matA(n,n) = eval_gaussian(parvec1,parvec1,parcovv);
  }
  double matAval = 0.0;
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      parvec1 = muA.row(i).t();
      parvec2 = muA.row(j).t();
      parcovv = covA.slice(i)+covA.slice(j);
      
      matAval = eval_gaussian(parvec1, parvec2, parcovv);
      matA(i,j) = matAval;
      matA(j,i) = matAval;
    }
  }
  
  // matB : M normals from mixture 2
  for (int m=0;m<M;m++){
    parvec2 = muB.row(m).t();
    parcovv = 2*covB.slice(m);
    matB(m,m) = eval_gaussian(parvec2,parvec2,parcovv);
  }
  double matBval = 0.0;
  for (int i=0;i<(M-1);i++){
    for (int j=(i+1);j<M;j++){
      parvec1 = muB.row(i).t();
      parvec2 = muB.row(j).t();
      parcovv = covB.slice(i) + covB.slice(j);
      
      matBval = eval_gaussian(parvec1, parvec2, parcovv);
      matB(i,j) = matBval;
      matB(j,i) = matBval;
    }
  }
  
  
  // matC : (N,M) cross of two mixtures
  double matCval = 0.0;
  for (int i=0;i<N;i++){
    parvec1 = muA.row(i).t();
    for (int j=0;j<M;j++){
      parvec2 = muB.row(j).t();
      parcovv = covA.slice(i) + covB.slice(j);
      matC(i,j) = eval_gaussian(parvec1, parvec2, parcovv);
    }
  }
  
  // return output
  return Rcpp::List::create(Rcpp::Named("A")=matA,
                            Rcpp::Named("B")=matB,
                            Rcpp::Named("C")=matC);
}

///////////////////////////////////////////////////////////////////
// 2. integrate_1d
// [[Rcpp::export]]
double integrate_1d(arma::vec &tseq, arma::vec &fval){
  // parameters
  int N = tseq.n_elem;
  double output = 0.0;
  double dt = 0.0;
  
  // iteration
  for (int i=0;i<(N-1);i++){
    dt = tseq(i+1)-tseq(i);
    output += (fval(i) + fval(i+1))*dt/2.0;
  }
  return(output);
}