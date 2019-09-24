#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//////////////////// Sub Routine
LogicalMatrix isweird(NumericMatrix x){
  const int n = x.nrow();
  LogicalMatrix out(n,n);
  
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      out(i,j) = ((x(i,j)==R_NegInf)||(x(i,j)==R_PosInf)||(NumericVector::is_na(x(i,j))));
    }
  }
  return out;
}
//////////////////// Main Routine
// [[Rcpp::export]]
Rcpp::NumericMatrix aux_shortestpath(NumericMatrix& wmat){
  // 3-1. get ready
  const int v = wmat.nrow();
  NumericMatrix dist(v,v);
  for (int i=0;i<v;i++){
    for (int j=0;j<v;j++){
      dist(i,j) = R_PosInf;
    }
  }
  // 3-2. initialization
  LogicalMatrix checker = isweird(wmat);
  
  // 3-3. Floyd-Warshall algorithm
  // 3-3-1. vertex
  for (int i=0;i<v;i++){
    dist(i,i) = 0;
  }
  // 3-3-2. edge list
  for (int i=0;i<v;i++){
    for (int j=0;j<v;j++){
      if (checker(i,j)==false){
        dist(i,j) = wmat(i,j);
      }
    }
  }
  // 3-3-3. main iteration
  for (int k=0;k<v;k++){
    for (int i=0;i<v;i++){
      for (int j=0;j<v;j++){
        if (dist(i,j)>(dist(i,k)+dist(k,j))){
          dist(i,j)=dist(i,k)+dist(k,j);
        }
      }
    }
  }
  
  // 3-4. return output
  return(dist);
}