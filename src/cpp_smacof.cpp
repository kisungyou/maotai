#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

/* Auxiliary Functions
 * init_by_cmds            : quick lazy implementation of CMDS
 * construct_Vinv_weighted : Vinv 
 * compute_raw_stress      : compute the raw stress
 * operation_B             : B(Z) operator
 */

arma::mat operation_B(arma::mat &DZ, arma::mat &D, arma::mat &W){
  int N = DZ.n_rows;
  arma::mat BZ(N,N,fill::zeros);
  
  // off-diagonals first
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      if (DZ(i,j)!=0){
        BZ(i,j) = -W(i,j)*D(i,j)/DZ(i,j);
        BZ(j,i) = BZ(i,j);
      }
    }
  }
  
  // diagoanls
  arma::rowvec rowBZ(N,fill::zeros);
  for (int i=0; i<N; i++){
    rowBZ = BZ.row(i);
    BZ(i,i) = -arma::accu(rowBZ);
  }
  
  return(BZ);
}

double compute_raw_stress(arma::mat &DZ, arma::mat &D, arma::mat &W){
  int N = DZ.n_rows;
  
  double output = 0.0;
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      output += std::pow(DZ(i,j) - D(i,j), 2.0)*W(i,j);
    }
  }
  return(output);
}


arma::mat init_by_cmds(arma::mat &D, int ndim){
  int N = D.n_rows;
  arma::mat D2 = arma::pow(D, 2.0);
  arma::mat J  = arma::eye<arma::mat>(N,N) - (arma::ones<arma::mat>(N,N)/(static_cast<double>(N)));
  arma::mat B  = -0.5*J*D2*J;  
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, B);
  
  arma::mat hey_mat = arma::fliplr(eigvec.tail_cols(ndim));
  arma::vec hey_vec = arma::sqrt(arma::reverse(eigval.tail(ndim)));
  
  arma::mat output = hey_mat*arma::diagmat(hey_vec);
  return(output);
}

arma::mat construct_Vinv_weighted(arma::mat &W){
  int N = W.n_rows;
  
  // compute V first
  arma::rowvec rowW(N,fill::zeros);
  arma::mat V(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      V(i,j) = -W(i,j);
      V(j,i) = V(i,j);
    }
  }
  
  for (int i=0; i<N; i++){
    rowW = W.row(i); rowW(i) = 0.0;
    V(i,i) = arma::accu(rowW);
  }
  
  // compute the inverse
  arma::mat Vinv = arma::inv(V + arma::eye<arma::mat>(N,N)) - ((arma::eye<arma::mat>(N,N))/static_cast<double>(N*N));
  return(Vinv);
}


// [[Rcpp::export]]
Rcpp::List src_smacof(arma::mat &D, arma::mat &W, int ndim, int maxiter, double abstol, bool use_gutman){
  // initialize via CMDS
  int N = D.n_rows;
  arma::mat old_X = init_by_cmds(D, ndim);
  arma::mat new_X(N, ndim, fill::zeros);
  
  arma::mat old_Xdist(N,N,fill::zeros);
  arma::mat new_Xdist(N,N,fill::zeros);
  for (int i=0; i<(N-1); i++){
    for (int j=(i+1); j<N; j++){
      old_Xdist(i,j) = arma::norm(old_X.row(i)-old_X.row(j), 2);
      old_Xdist(j,i) = old_Xdist(i,j);
    }
  }
  
  // prepare
  double old_cost = compute_raw_stress(old_Xdist, D, W);
  double new_cost = 0.0;
  double inc_cost = 10000.0;
  
  arma::mat BZ(N,N,fill::zeros);
  arma::mat Vinv(N,N,fill::zeros);
  if (use_gutman){
    Vinv = arma::eye<arma::mat>(N,N)/(static_cast<double>(N));
  } else {
    Vinv = construct_Vinv_weighted(W);
  }
      
  // iterate
  for (int it=0; it<maxiter; it++){
    // compute the iterate
    BZ = operation_B(old_Xdist, D, W);
    if (use_gutman){
      new_X = Vinv*BZ*old_X;
    } else {
      new_X = BZ*old_X/(static_cast<double>(N));
    }
    // compute the pairwise distance
    for (int i=0; i<(N-1); i++){
      for (int j=(i+1); j<N; j++){
        new_Xdist(i,j) = arma::norm(new_X.row(i)-new_X.row(j), 2);
        new_Xdist(j,i) = new_Xdist(i,j);
      }
    }
    new_cost = compute_raw_stress(new_Xdist, D, W);
    inc_cost = std::abs(new_cost - old_cost);
    
    old_cost  = new_cost;
    old_X     = new_X;
    old_Xdist = new_Xdist;
    if (inc_cost < abstol){
      break;
    }
  }
  
  
  return Rcpp::List::create(Rcpp::Named("embed")=old_X, 
                            Rcpp::Named("stress")=old_cost);
}
