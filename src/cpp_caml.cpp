#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// -------------------------------------------------------------------------
// HELPER FUNCTIONS
// calm_skew
//    - skew(A) = (A - A^T)/2
//    - extract vertical component
// calm_lyapunov_skew(S, K)
//    - solve SX + XS = K
//    - equivalent to SX + XS - K = 0
// caml_qf_qr
//    - thin QR factorization Y=QR and returns the Q factor
// caml_procrustes_align
//    - Find optimal Q such that |Y-Yi*Q| is minimized
// caml_proj_horizontal
//    - Proj_Y (Z) = Z - YX, SX+XS = skew(Y^T Z)
// caml_dist2_quotient
//    - squared quotient distance
// -------------------------------------------------------------------------
arma::mat caml_skew(const arma::mat& A){
  return(0.5*(A - A.t()));
}
arma::mat caml_lyapunov_skew(const arma::mat& S, const arma::mat& K) {
  arma::mat Omega;
  
  // Armadillo solves: A*X + X*B + C = 0
  // We want: S*Omega + Omega*S = K
  // => S*Omega + Omega*S - K = 0, so C = -K
  bool ok = arma::sylvester(Omega, S, S, -K);
  
  if(!ok) {
    Rcpp::stop("caml_lyapunov_skew: sylvester solver failed (S may be ill-conditioned).");
  }
  
  // enforce skew-symmetry numerically
  Omega = caml_skew(Omega);
  return Omega;
}
arma::mat caml_qf_qr(const arma::mat& Y) {
  arma::mat Q, R;
  arma::qr_econ(Q, R, Y); // thin QR
  
  // Fix sign ambiguity so diag(R) is positive
  arma::vec d = arma::diagvec(R);
  arma::vec s = arma::sign(d);
  s.elem(arma::find(s == 0)).ones();
  Q = Q * arma::diagmat(s);
  
  return Q; // n x k with orthonormal columns
}
arma::mat caml_procrustes_align(const arma::mat& Y, const arma::mat& Yi) {
  arma::mat A = Yi.t() * Y;   // k x k
  arma::mat U, V;
  arma::vec s;
  arma::svd(U, s, V, A);      // A = U diag(s) V^T
  arma::mat Q = U * V.t();    // optimal orthogonal
  return Yi * Q;              // aligned factor
}
arma::mat caml_proj_horizontal(const arma::mat& Y, const arma::mat& Z) {
  arma::mat S = Y.t() * Y;            // k x k SPD
  arma::mat K = caml_skew(Y.t() * Z); // k x k skew
  arma::mat Omega = caml_lyapunov_skew(S, K);
  return Z - Y * Omega;              // horizontal projection
}
double caml_dist2_quotient(const arma::mat& Y, const arma::mat& Yi) {
  arma::mat Yi_al = caml_procrustes_align(Y, Yi);
  arma::mat D = Y - Yi_al;
  return arma::accu(D % D); // ||D||_F^2
}
arma::mat caml_centering(const arma::mat& X){
  int N = X.n_rows;
  int P = X.n_cols;
  arma::rowvec Xmean = arma::mean(X, 0);
  arma::mat output(N,P,fill::zeros);
  for (int n=0; n<N; n++){
    output.row(n) = X.row(n) - Xmean;
  }
  return(output);
}
arma::mat caml_retr_gauge(const arma::mat& Y, const arma::mat& eta) {
  arma::mat Ytmp = Y + eta;
  
  arma::mat Q, R;
  arma::qr_econ(Q, R, Ytmp);
  
  arma::vec d = arma::diagvec(R);
  arma::vec s = arma::sign(d);
  s.elem(arma::find(s == 0)).ones();          // avoid zeros
  arma::mat Qsign = arma::diagmat(s);         // orthogonal (±1)
  
  // Right-multiply by orthogonal Qsign: preserves Ytmp Ytmp^T
  return Ytmp * Qsign;
}
arma::mat caml_gauge_sign(const arma::mat& Y){
  arma::mat Q, R;
  arma::qr_econ(Q, R, Y);
  
  arma::vec d = arma::diagvec(R);
  arma::vec s = arma::sign(d);
  s.elem(arma::find(s == 0)).ones();
  arma::mat Qsign = arma::diagmat(s);
  
  return Y * Qsign;  // preserves YY^T
}



// FRECHET MEAN COMPUTATION ---------------------------------------------------
double caml_frechet_cost(const arma::mat& Y,
                         const std::vector<arma::mat>& Ylist,
                         const arma::vec& w) {
  double val = 0.0;
  for (size_t i = 0; i < Ylist.size(); i++) {
    val += w(i) * caml_dist2_quotient(Y, Ylist[i]);
  }
  return 0.5 * val;
}
arma::mat caml_frechet_egrad(const arma::mat& Y,
                             const std::vector<arma::mat>& Ylist,
                             const arma::vec& w) {
  int n = Y.n_rows;
  int p = Y.n_cols;
  
  arma::mat G(n,p,fill::zeros);
  for (size_t i = 0; i < Ylist.size(); i++) {
    arma::mat Yi_al = caml_procrustes_align(Y, Ylist[i]);
    G += w(i) * (Y - Yi_al);
  }
  return G;
}

// [[Rcpp::export]]
Rcpp::List caml_frechet_mean(const Rcpp::List& Y_list_R,
                            const arma::vec& w,
                            int maxit = 200,
                            double tol = 1e-8,
                            double step0 = 1.0,
                            double c1 = 1e-4,
                            double beta = 0.5,
                            bool verbose = false) {
  
  int N = Y_list_R.size();
  std::vector<arma::mat> Ylist(N);
  arma::mat Y_tmp;
  for(int i = 0; i < N; i++){ // apply centering
    Y_tmp.reset();
    Y_tmp = Rcpp::as<arma::mat>(Y_list_R[i]);
    Ylist[i] = caml_centering(Y_tmp);
  } 
  
  if((int)w.n_elem != N) Rcpp::stop("w must have length N");
  if(std::abs(arma::accu(w) - 1.0) > 1e-6) Rcpp::stop("weights w must sum to 1");
  
  // init: align to first, average, then qf
  arma::mat Y = w(0) * Ylist[0];
  for(int i=1; i<N; i++){
    arma::mat Yi_al = caml_procrustes_align(Ylist[0], Ylist[i]);
    Y += w(i) * Yi_al;
  }
  Y = caml_gauge_sign(Y);
  
  double f0 = caml_frechet_cost(Y, Ylist, w);
  if(verbose) Rcpp::Rcout << "init cost = " << f0 << "\n";
  
  for(int it = 0; it < maxit; it++) {
    
    arma::mat eG = caml_frechet_egrad(Y, Ylist, w);
    arma::mat rG = caml_proj_horizontal(Y, eG);
    double gnorm = arma::norm(rG, "fro");
    
    if(verbose) Rcpp::Rcout << "it " << it
                            << " cost " << f0
                            << " |grad| " << gnorm << "\n";
    
    if(gnorm < tol) break;
    
    arma::mat eta = -rG;  // descent direction
    
    // Armijo backtracking along retraction curve
    double step = step0, fnew;
    arma::mat Ynew;
    while(true) {
      // Ynew = caml_qf_qr(Y + step * eta);
      Ynew = caml_retr_gauge(Y, step * eta);
      fnew = caml_frechet_cost(Ynew, Ylist, w);
      if(fnew <= f0 - c1 * step * gnorm * gnorm) break;
      step *= beta;
      if(step < 1e-16) break;
    }
    
    if(std::abs(f0 - fnew) < tol * (1.0 + f0)) {
      Y = Ynew; f0 = fnew; break;
    }
    Y = Ynew; f0 = fnew;
  }
  
  arma::mat X = Y * Y.t();
  return Rcpp::List::create(
    Rcpp::_["Y"] = Y,
    Rcpp::_["YYT"] = X,
    Rcpp::_["cost"] = f0
  );
}
