#include <RcppDist.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

double my_invgamma(double alpha, double beta){
  return(1.0/R::rgamma(alpha,1.0/beta));
}

// 2. my_dinvgamma : inverse gamma evaluator
double my_dinvgamma(double x, double alpha, double beta){
  return(1.0/R::dgamma(x, alpha, 1.0/beta, 0));
}

// Auxiliary Function : compute SSR
// [[Rcpp::export]]
double compute_SSR(arma::mat &D, arma::mat &Delta){
  // parameters
  int N = D.n_rows;
  double NN = static_cast<double>(N);
  
  // compute via iteration
  double outval = 0.0;
  double tobesq = 0.0;
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tobesq = (D(i,j)-Delta(i,j));
      outval += (tobesq*tobesq)/NN;
    }
  }
  return(outval);
}
double compute_SSR_xmat(arma::mat &D, arma::mat &Xnew){ // this one is using matrix data
  int N = D.n_rows; double NN = static_cast<double>(N);
  int p = Xnew.n_cols;
  
  double outval = 0.0;
  double tobesq = 0.0;
  
  arma::rowvec xvec1(p,fill::zeros);
  arma::rowvec xvec2(p,fill::zeros);
  
  double Delij = 0.0;
  for (int i=0;i<N;i++){
    xvec1 = Xnew.row(i);
    for (int j=(i+1);j<N;j++){
      xvec2  = Xnew.row(j);
      Delij  = arma::norm(xvec1-xvec2, 2);
      tobesq = D(i,j)-Delij;
      outval+= (tobesq*tobesq)/NN;
    }
  }
  return(outval);
}
// Auxiliary Function : compute pairwise distance function
arma::mat compute_pdmat(arma::mat &X){
  int N = X.n_rows;
  int p = X.n_cols;
  arma::mat output(N,N,fill::zeros);
  arma::vec tgt1(p,fill::zeros);
  arma::vec tgt2(p,fill::zeros);
  double tmpval = 0.0;
  for (int i=0;i<(N-1);i++){
    tgt1 = X.row(i).t();
    for (int j=0;j<N;j++){
      tgt2 = X.row(j).t();
      tmpval = arma::norm(tgt1-tgt2,2);
      output(i,j) = tmpval;
      output(j,i) = tmpval;
    }
  }
  return(output);
}
// Auxiliary Function : centering and rotating
arma::mat crotX(arma::mat &X){
  int N = X.n_rows;
  int p = X.n_cols;
  
  arma::mat Xtmp(N,p,fill::zeros);
  arma::rowvec xmean = arma::mean(X, 0);
  for (int i=0;i<N;i++){
    Xtmp.row(i) = X.row(i)-xmean;
  }
  
  arma::mat Xcov = Xtmp.t()*Xtmp/(static_cast<double>(N));
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Xcov);
  
  arma::mat output = Xtmp*eigvec;
  return(output);
}
// Auxiliary Function : a single step for MH update
arma::rowvec update_xvec(arma::mat D, arma::mat X, int id, double sigma2, double constant, arma::mat Lbdmat){
  int N = X.n_rows; double NN = static_cast<double>(N);
  int p = X.n_cols;
  arma::mat Xold = X;
  arma::mat Xtgt = X;
  double stepsize = static_cast<double>(std::sqrt(static_cast<float>(sigma2*constant/(NN-1.0))));
  for (int i=0;i<p;i++){
    Xtgt(id,i) += R::rnorm(0.0, stepsize);
  }
  double sigma = sqrt(sigma2);
  
  arma::vec xtgt = Xtgt.row(id).t(); // column vectors
  arma::vec xold = Xold.row(id).t();
  
  // common variables
  double tmpval = 0.0;
  
  // need to evaluate two ratio
  // (1) compute for xtgt
  arma::mat Deltgt = compute_pdmat(Xtgt);
  double Q1tgt = 0.0; 
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tmpval = D(i,j)-Deltgt(i,j);
      Q1tgt += (tmpval*tmpval)/sigma2;
    }
  }
  double Q2tgt = arma::dot(xtgt, arma::solve(Lbdmat, xtgt)); 
  double t3tgt = 0.0;
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      if (i!=j){
        t3tgt += static_cast<double>(std::sqrt(static_cast<float>(R::pnorm5(Deltgt(i,j)/sigma,0.0,1.0,1,0))));
      }
    }
  }
  double ftgt = -(Q1tgt+Q2tgt)/2.0 - t3tgt;
  
  // (2) compute for xold
  arma::mat Delold = compute_pdmat(Xold);
  double Q1old = 0.0; 
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tmpval = D(i,j)-Delold(i,j);
      Q1old += (tmpval*tmpval)/sigma2;
    }
  }
  double Q2old = arma::dot(xold, arma::solve(Lbdmat, xold)); 
  double t3old = 0.0;
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      if (i!=j){
        t3old += static_cast<double>(std::sqrt(static_cast<float>(R::pnorm5(Delold(i,j)/sigma,0.0,1.0,1,0))));
      }
    }
  }
  double fold = -(Q1old+Q2old)/2.0 - t3old;
  
  // (3) compute the ratio (?)
  double fratio = exp(ftgt-fold);
  if (fratio >= 1){
    fratio = 1.0;
  }
  double rnumbr = R::runif(0.0, 1.0);
  if (rnumbr <= fratio){ // accept
    return(xtgt.t());
  } else {
    return(xold.t());
  }
} 

// Auxiliary Function : 'stress'
// https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Multidimensional_Scaling.pdf
// [[Rcpp::export]]
double compute_stress(arma::mat &D, arma::mat &Dhat){ // D is original distance, Dhat is estimated ones
  int N = D.n_rows;
  
  double tobesq = 0.0;
  double term1  = 0.0; // numerator
  double term2  = 0.0; // denominator
  for (int i=0;i<(N-1);i++){
    for (int j=(i+1);j<N;j++){
      tobesq = D(i,j)-Dhat(i,j);
      term1 += (tobesq*tobesq);
      term2 += D(i,j)*D(i,j);
    }
  }
  return(sqrt(term1/term2));
}

// Main Computation
// D  : given pairwise distances
// X0 : (initial) MDS locations
// 
// [[Rcpp::export]]
Rcpp::List main_bmds(arma::mat D, arma::mat X0, double sigg0, 
                     double a, double alpha, int maxiter, double constant, bool verbose, 
                     arma::vec betas){
  // 1) some parameters
  int N = X0.n_rows; double NN = static_cast<double>(N);
  int p = X0.n_cols;
  double m = NN*(NN-1.0)/2.0;
  
  
  // 2) setup 
  arma::mat Xold = crotX(X0); // X will not be recorded, just use 
  arma::mat Xnew(N,p,fill::zeros); 
  arma::mat Xsol = Xold;
  
  double SSRnew = 0.0;
  double SSRold = compute_SSR_xmat(D, Xold);
  double SSRsol = SSRold;
  
  arma::mat Sold(p,p,fill::zeros);
  double sigma2 = sigg0;
  double sigtmp = 0.0;
  arma::vec vecs(p,fill::zeros);
  arma::vec lambdas(p,fill::zeros);
  arma::mat Lbdmat;
  arma::rowvec tmprow(p,fill::zeros);
  double b = (a-1)*SSRold/m; // paper's setup
  
  double varalpha = 0.0;
  double varbeta  = 0.0;
  double varvar = 0.0;
  double varratio = 0.0;
  
  // 3) iteration
  // int accept = 0;
  for (int i=0;i<maxiter;i++){
    // 3-1. update lambdas
    for (int j=0;j<p;j++){ // compute sample variances for each coordinate
      vecs(j) = arma::var(Xold.col(j))*NN;
    }
    Sold = Xold.t()*Xold/NN;
    for (int j=0;j<p;j++){ // sample from IG
      lambdas(j) = my_invgamma(alpha+NN/2.0, betas(j) + vecs(j)/2.0); // according to the paper's choice
    }
    Lbdmat = arma::diagmat(lambdas);
    
    // 3-2. update X 
    Xnew   = Xold;
    for (int j=0;j<N;j++){ // for each row
      tmprow = update_xvec(D, Xnew, j, sigma2, constant, Lbdmat);
      Xnew.row(j) = tmprow;
    }
    SSRnew = compute_SSR_xmat(D, Xnew); // update SSR
    Xnew   = crotX(Xnew); // centering + rotation
    
    // 3-3. update sigma using MH
    varalpha = m/2 + a;
    varbeta  = SSRnew/2 + b;
    varvar   = (varbeta*varbeta)/((varalpha-1)*(varalpha-1)*(varalpha-2));
    
    sigtmp = sigma2 + R::rnorm(0, sqrt(constant*varvar));
    if (sigtmp > 0){ // let's compare
      varratio = my_dinvgamma(sigtmp,varalpha,varbeta)/my_dinvgamma(sigma2,varalpha,varbeta);
      if (varratio > 1){
        varratio = 1.0;
      }
      if (R::runif(0,1) <= varratio){
        sigma2 = sigtmp;
      }
    }
    
    
    // 3-4. update correspondingly
    if (SSRnew < SSRsol){ // running record of the best solution
      SSRsol = SSRnew;
      Xsol   = Xnew;
    }
    SSRold = SSRnew;
    Xold   = Xnew;
    
    // 3-5. report the update
    if (verbose==true){
      Rcpp::Rcout << "** bmds : iteration " << i+1 << "/" << maxiter << " complete." << std::endl;
    }
  }
  
  // 4) return
  return Rcpp::List::create(Rcpp::Named("solX")=Xsol);
}
