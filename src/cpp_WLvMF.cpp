#include <RcppArmadillo.h>
#include "math.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// basic operations
arma::rowvec sphere_exp(arma::rowvec x, arma::rowvec d, double t){
  double nrm_td = arma::norm(t*d, 2);
  int N = x.n_elem;
  arma::rowvec out(N,fill::zeros);
  if (nrm_td < 1e-15){
    out = x;
  } else {
    arma::rowvec tmp(N,fill::zeros);
    tmp = std::cos(nrm_td)*x + ((sin(nrm_td))/nrm_td)*t*d;
    out = tmp/arma::norm(tmp, 2);
  }
  return(out);
}

arma::rowvec sphere_proj(arma::rowvec x, arma::rowvec u){
  return(u-x*(arma::dot(x,u)));
}
// [[Rcpp::export]]
double sphere_dist(arma::rowvec vecx, arma::rowvec vecy){
  arma::rowvec vecxy = vecx-vecy;
  double dotxy = arma::dot(vecx, vecy);
  
  if (arma::norm(vecxy, 2) < arma::datum::eps){
    return(0.0);
  } else if (std::sqrt(dotxy*dotxy) >= (1.0-arma::datum::eps)){
    return(arma::datum::pi);
  } else {
    return(std::acos(arma::dot(vecx, vecy)));  
  }
}
arma::rowvec sphere_log(arma::rowvec x, arma::rowvec y){
  arma::rowvec v = sphere_proj(x,y-x);
  double di = sphere_dist(x,y);
  if (di > 1e-6){
    double nv = arma::norm(v, 2);
    v = v*(di/nv);
  }
  return(v);
}

// MAIN ROUTINES ===============================================================
// row-normalization

// [[Rcpp::export]]
arma::mat cpp_WL_normalise(arma::mat& X){
  int N = X.n_rows;
  int d = X.n_cols;
  
  arma::mat output(N,d,fill::zeros);
  for (int n=0; n<N; n++){
    output.row(n) = X.row(n)/arma::norm(X.row(n),2);
  }
  return(output);
}

// weighted mean computation
// // [[Rcpp::export]]
// arma::rowvec cpp_WL_weighted_mean(arma::mat& X, arma::vec& weights){
//   int N = X.n_rows;
//   int d = X.n_cols;
//   
//   arma::rowvec Sold = arma::mean(X, 0);
//   double Sinc = 0.0;
//   Sold = Sold/arma::norm(Sold,2);
//   arma::rowvec distvec(N,fill::zeros);
//   arma::rowvec Snew(d,fill::zeros);
//   arma::rowvec Stmp(d,fill::zeros);
//   
//   for (int it=0; it<1000; it++){
//     // reset the temporary gradient 
//     Stmp.fill(0.0);
//     for (int n=0; n<N; n++){
//       Stmp += 2.0*weights(n)*sphere_log(Sold, X.row(n));
//     }
//     Snew = sphere_exp(Sold, Stmp, 1.0);
//     Sinc = arma::norm(Sold-Snew,2);
//     Sold = Snew;
//     
//     if (Sinc < 1e-7){
//       break;
//     }
//   }
//   
//   return(Sold);
// }
// [[Rcpp::export]]
arma::rowvec cpp_WL_weighted_mean(arma::mat& X, arma::vec& weights){
  int N = X.n_rows;
  int d = X.n_cols;
  
  arma::rowvec Sout(d, fill::zeros);
  for (int n=0; n<N; n++){
    Sout += (weights(n)/arma::accu(weights))*X.row(n);
  }
  
  arma::rowvec Sfin = Sout/arma::norm(Sout, 2);
  return(Sfin);
}
