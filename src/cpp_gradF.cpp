#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat gradF(Function func, arma::mat xnow, double h){
  int m = xnow.n_rows;
  int n = xnow.n_cols;
  arma::mat dX(m,n,fill::zeros);
  arma::mat Xp = xnow;
  arma::mat Xm = xnow;
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      Xp(i,j) = Xp(i,j) + h;
      Xm(i,j) = Xm(i,j) - h;
      dX(i,j) = (sum(as<NumericVector>(func(Xp)))- sum(as<NumericVector>(func(Xm))))/(2.0*h);
      Xp(i,j) = Xp(i,j) - h;
      Xm(i,j) = Xm(i,j) + h;
    }
  }
  
  return(dX);
}