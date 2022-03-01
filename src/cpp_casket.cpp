#include <RcppArmadillo.h>
#include "evaluations.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * 1. aux_shortestpath  : 'shortestpath'
 * 2. cppsub_2007Wang   : 'trio'
 * 3. gradF             : 'matderiv'
 * 4. dat2centers       : 'dpmeans'
 * 5. cpp_sylvester &
 *    cpp_lyapunov
 * 6. cpp_weiszfeld     : 'weiszfeld'
 * 7. cpp_kmeans        : 'kmeans' for testing only
 * 8. emds_gamma0       : 'emds'
 * 9. cpp_pairwise_L2 : L2 distance between GMM's.
 * 10. integrate_1d    : 1d integration, generic but used in distance computation
 * 11. cpp_pdist       : compute pairwise distance
 * 12. cpp_geigen      : do 'geigen' pairwise eigendecomposition
 * 13. cpp_triangle    : check triangle inequality
 */

///////////////////////////////////////////////////////////////////
// 1. aux_shortestpath
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

///////////////////////////////////////////////////////////////////
// 2. cppsub_2007Wang
// [[Rcpp::export]]
arma::mat cppsub_2007Wang(arma::mat V0, int mm, int d, arma::mat Spu, arma::mat Stu, int maxiter, double eps){
  // 1. preliminary setup
  double abstol = std::sqrt((static_cast<double>(mm*d))*eps);
  
  arma::mat Vold = V0;
  arma::mat Vnew(mm,d,fill::zeros);
  
  double lbdn = 0.0;
  double incV = 0.0;
  
  arma::vec Vval(mm,fill::zeros);
  arma::mat Vvec(mm,mm,fill::zeros);
  arma::mat Vtmp(mm,d,fill::zeros);
  
  arma::mat Stv(mm,mm,fill::zeros);
  
  // 2. do the iteration
  for (int i=0;i<maxiter;i++){
    // 2-1. compute the lambda
    lbdn = arma::trace(Vold.t()*Spu*Vold)/arma::trace(Vold.t()*Stu*Vold);
    // 2-2. solve the eigenvalue problem
    eig_sym(Vval,Vvec,(Spu-(lbdn*Stu)));
    // 2-3. extract Vtmp
    Vtmp = Vvec.tail_cols(d);
    // 2-4. update by readjusting
    Stv  = Vtmp*Vtmp.t()*Stu*Vtmp*Vtmp.t();
    eig_sym(Vval,Vvec,Stv);
    Vnew = Vvec.tail_cols(d);
    // 2-5. updating info and update
    incV = arma::norm(Vold-Vnew,"fro");
    Vold = Vnew;
    // 2-6. stop if the criterion is met
    if (incV < abstol){
      break;
    }
    // no need to print; Rcpp::Rcout << "* iteration " << i+1 << " complete.." << std::endl;
  }
  
  // 3. return the value
  return(Vold);
}

///////////////////////////////////////////////////////////////////
// 3. gradF
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

///////////////////////////////////////////////////////////////////
// 4. dat2centers
// [[Rcpp::export]]
arma::vec dat2centers(arma::rowvec data, arma::mat &centers){
  // parameters
  int K = centers.n_rows;
  int p = data.n_cols;
  
  // compute
  arma::vec dic(K,fill::zeros);
  arma::rowvec diffvec(p,fill::zeros);
  for (int k=0;k<K;k++){
    diffvec = data-centers.row(k);
    dic(k)  = arma::dot(diffvec,diffvec);
  }
  
  // report
  return(dic);
}

///////////////////////////////////////////////////////////////////
// 5. cpp_sylvester & cpp_lyapunov
// [[Rcpp::export]]
arma::mat cpp_sylvester(arma::mat A, arma::mat B, arma::mat C){
  arma::mat solution;
  arma::syl(solution,A,B,C);
  return(solution);
}

// [[Rcpp::export]]
arma::mat solve_lyapunov(arma::mat A, arma::mat B, arma::mat C){
  // simply solve it !
  arma::mat solution;
  arma::syl(solution, A, B, C);
  return(solution);
}

///////////////////////////////////////////////////////////////////
// 6. weiszfeld
// [[Rcpp::export]]
arma::rowvec cpp_weiszfeld(arma::mat X, double abstol, int maxiter, arma::rowvec xinit, arma::vec weights, double epsnum){
  // parameters
  int N = X.n_rows;
  int P = X.n_cols;
  
  // prepare
  arma::rowvec xold = xinit;
  arma::rowvec xtmp(P,fill::zeros);
  arma::rowvec xnew(P,fill::zeros);
  arma::vec dists(N,fill::zeros);
  double xtmp2 = 0.0;
  double xinc  = 0.0;
  double norm2 = 0.0;
  
  // iteration
  for (int it=0;it<maxiter;it++){
    // step 1. compute distance
    for (int n=0;n<N;n++){
      norm2 = arma::norm(X.row(n)-xold, 2);
      if (norm2 < epsnum){
        dists(n) = norm2 + epsnum;
      } else {
        dists(n) = norm2;
      }
    }
    // step 2. compute numerator and denominator
    xtmp.fill(0.0);
    xtmp2 = 0.0;
    for (int n=0;n<N;n++){
      xtmp  += weights(n)*X.row(n)/dists(n);
      xtmp2 += weights(n)/dists(n);
    }
    xnew = xtmp/xtmp2;
    
    // step 3. updating information
    xinc = arma::norm(xold-xnew,2);
    xold = xnew;
    if (xinc < abstol){
      break;
    }
  }
  
  // return
  return(xold);
}

///////////////////////////////////////////////////////////////////
// 7. kmeans test from RcppArmadillo implementation
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

///////////////////////////////////////////////////////////////////
// 8. compute gamma0 for EMDS
// [[Rcpp::export]]
double emds_gamma0(arma::mat dmat){
  // parameters
  int N = dmat.n_rows;
  double gamma0 = 0.0;
  double theval = 0.0;
  
  // iterate.. triplet
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      for (int k=0;k<N;k++){
        theval = std::abs(dmat(i,j)+dmat(i,k)-dmat(j,k));
        if (theval > gamma0){
          gamma0 = theval;
        }
      }
    }
  }
  
  // report
  return(gamma0);
}

///////////////////////////////////////////////////////////////////
// 9. cpp_pairwise_L2
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
  // double matCval = 0.0;
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
// 10. integrate_1d
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
///////////////////////////////////////////////////////////////////
// 11. cpp_pdist
// [[Rcpp::export]]
arma::mat cpp_pdist(arma::mat X){
  // parameters
  int N = X.n_rows;
  int d = X.n_cols;
  
  // prepare
  arma::mat output(N,N,fill::zeros);
  arma::rowvec xdiff(d,fill::zeros);
  // iteration
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      xdiff = X.row(i)-X.row(j);
      output(i,j) = arma::norm(xdiff, 2);
      output(j,i) = output(i,j);
    }
  }
  return(output);
}

///////////////////////////////////////////////////////////////////
// 12. cpp_geigen
// [[Rcpp::export]]
Rcpp::List cpp_geigen(arma::mat& A, arma::mat& B){
  arma::cx_vec eigval;
  arma::cx_mat eigmat;
  
  eig_pair(eigval, eigmat, A, B);
  return Rcpp::List::create(Rcpp::Named("values")=eigval,
                            Rcpp::Named("vectors")=eigmat);
}

///////////////////////////////////////////////////////////////////
// 13. cpp_triangle
// [[Rcpp::export]]
bool cpp_triangle(arma::mat& D){
  int N = D.n_rows;
  
  double term1 = 0.0;
  double term2 = 0.0;
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      term1 = D(i,j);
      for (int k=0;k<N;k++){
        term2 = D(i,k) + D(k,j) + arma::datum::eps;
        if (!(term1 <= term2)){
          return(false);
        }
      }
    }
  }
  return(true);
}
