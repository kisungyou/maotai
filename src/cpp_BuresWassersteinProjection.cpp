// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;

// --------- small helpers ---------

static inline mat symm(const mat& A) {
  return 0.5 * (A + A.t());
}

static inline mat spd_sqrt(const mat& A, vec* evals_out = nullptr, mat* evecs_out = nullptr, double eps = 1e-12) {
  vec evals;
  mat evecs;
  arma::eig_sym(evals, evecs, symm(A));
  evals = arma::clamp(evals, eps, arma::datum::inf);
  if (evals_out) *evals_out = evals;
  if (evecs_out) *evecs_out = evecs;
  vec s = arma::sqrt(evals);
  return evecs * arma::diagmat(s) * evecs.t();
}

// Solve Sylvester: S X + X S = M for symmetric SPD S (self-adjoint operator)
static inline mat sylv_self_adjoint_solve(const mat& S, const mat& M, double eps = 1e-12) {
  vec evals;
  mat V;
  arma::eig_sym(evals, V, symm(S));
  evals = arma::clamp(evals, eps, arma::datum::inf);
  
  mat Mt = V.t() * M * V;
  const uword d = evals.n_elem;
  mat Xt(d, d, fill::zeros);
  
  for (uword a = 0; a < d; ++a) {
    for (uword b = 0; b < d; ++b) {
      Xt(a,b) = Mt(a,b) / (evals(a) + evals(b));
    }
  }
  mat X = V * Xt * V.t();
  return symm(X);
}

// Gaussian W2^2 between N(m1,S1) and N(m2,S2)
static inline double gauss_w2_sq(const vec& m1, const mat& S1,
                                 const vec& m2, const mat& S2,
                                 double eps = 1e-12) {
  vec dm = m1 - m2;
  double mean_term = dot(dm, dm);
  
  mat S2_sqrt = spd_sqrt(S2, nullptr, nullptr, eps);
  mat C = S2_sqrt * S1 * S2_sqrt;
  mat C_sqrt = spd_sqrt(C, nullptr, nullptr, eps);
  
  double bures = trace(S1) + trace(S2) - 2.0 * trace(C_sqrt);
  return mean_term + bures;
}

// Projected W2^2 and its Euclidean gradient w.r.t U for a pair (i,j).
static inline void pair_s_and_gradU(
    const mat& U,
    const vec& mi, const mat& Si,
    const vec& mj, const mat& Sj,
    double& s_out,
    mat& gradU_out,
    double eps = 1e-12
) {
  const uword p = U.n_rows;
  const uword d = U.n_cols;
  
  vec a = mi - mj;
  
  // Mean part: || U^T a ||^2
  vec Ut_a = U.t() * a;
  double mean_part = dot(Ut_a, Ut_a);
  
  // Covariance part: h(A,B)
  mat A = symm(U.t() * Si * U);
  mat B = symm(U.t() * Sj * U);
  A.diag() += eps;
  B.diag() += eps;
  
  mat B_sqrt = spd_sqrt(B, nullptr, nullptr, eps);
  mat C = symm(B_sqrt * A * B_sqrt);
  mat S = spd_sqrt(C, nullptr, nullptr, eps); // S = C^{1/2}
  
  mat I = eye<mat>(d,d);
  mat G = sylv_self_adjoint_solve(S, I, eps);
  
  mat gradA_h = I - 2.0 * symm(B_sqrt * G * B_sqrt);
  
  mat M = A * G * B_sqrt + B_sqrt * G * A;
  mat Z = sylv_self_adjoint_solve(B_sqrt, M, eps);
  mat gradB_h = I - 2.0 * symm(Z);
  
  double cov_part = trace(A) + trace(B) - 2.0 * trace(S);
  
  s_out = mean_part + cov_part;
  
  gradU_out.set_size(p, d);
  gradU_out = 2.0 * (a * a.t()) * U
  + 2.0 * Si * U * gradA_h
  + 2.0 * Sj * U * gradB_h;
}

// Same projected W2^2 but WITHOUT gradient (for line search objective eval)
static inline double pair_s_only(
    const mat& U,
    const vec& mi, const mat& Si,
    const vec& mj, const mat& Sj,
    double eps = 1e-12
) {
  vec a = mi - mj;
  vec Ut_a = U.t() * a;
  double mean_part = dot(Ut_a, Ut_a);
  
  mat A = symm(U.t() * Si * U);
  mat B = symm(U.t() * Sj * U);
  A.diag() += eps;
  B.diag() += eps;
  
  mat B_sqrt = spd_sqrt(B, nullptr, nullptr, eps);
  mat C = symm(B_sqrt * A * B_sqrt);
  mat S = spd_sqrt(C, nullptr, nullptr, eps);
  
  double cov_part = trace(A) + trace(B) - 2.0 * trace(S);
  return mean_part + cov_part;
}

// QR retraction
static inline mat qr_retraction(const mat& X) {
  mat Q, R;
  qr_econ(Q, R, X);
  vec sgn = sign(R.diag());
  for (uword k = 0; k < sgn.n_elem; ++k) if (sgn(k) == 0) sgn(k) = 1.0;
  Q = Q * diagmat(sgn);
  return Q;
}

// Objective only (weights=1, full pairs)
static inline double objective_only(
    const mat& U,
    const mat& M,
    const cube& Sigma,
    const mat& D,
    double eps = 1e-12
) {
  const uword n = M.n_rows;
  double obj = 0.0;
  
  for (uword i = 0; i < n; ++i) {
    vec mi = M.row(i).t();
    mat Si = symm(Sigma.slice(i));
    Si.diag() += eps;
    
    for (uword j = i+1; j < n; ++j) {
      vec mj = M.row(j).t();
      mat Sj = symm(Sigma.slice(j));
      Sj.diag() += eps;
      
      double s_ij = pair_s_only(U, mi, Si, mj, Sj, eps);
      double delta = std::sqrt(std::max(s_ij, 0.0) + eps);
      double r = delta - D(i,j);
      obj += r * r;
    }
  }
  return obj;
}

// [[Rcpp::export]]
Rcpp::List bwp_project_gaussians_rcpp(
    const arma::mat& M,          // n x p
    const arma::cube& Sigma,     // p x p x n
    const int d,
    const int max_iter = 200,
    const double step_init = 1.0,     // initial step for Armijo
    const double tol_grad = 1e-6,
    const double tol_obj  = 1e-10,
    const double eps = 1e-12,
    const bool verbose = false,
    const double armijo_c = 1e-4,     // sufficient decrease constant
    const double armijo_beta = 0.5,   // backtracking shrink factor in (0,1)
    const int ls_max = 25,            // max backtracking steps
    const double ls_min_step = 1e-12  // minimum allowed step
) {
  const uword n = M.n_rows;
  const uword p = M.n_cols;
  
  if ((uword)d > p) Rcpp::stop("d must be <= p.");
  if (Sigma.n_rows != p || Sigma.n_cols != p || Sigma.n_slices != n) {
    Rcpp::stop("Sigma must have size p x p x n (R dim=c(p,p,n)).");
  }
  if (step_init <= 0) Rcpp::stop("step_init must be positive.");
  if (armijo_beta <= 0 || armijo_beta >= 1) Rcpp::stop("armijo_beta must be in (0,1).");
  if (armijo_c <= 0 || armijo_c >= 1) Rcpp::stop("armijo_c must be in (0,1).");
  
  // Precompute high-D pairwise distances D
  mat D(n, n, fill::zeros);
  for (uword i = 0; i < n; ++i) {
    vec mi = M.row(i).t();
    mat Si = symm(Sigma.slice(i));
    Si.diag() += eps;
    
    for (uword j = i+1; j < n; ++j) {
      vec mj = M.row(j).t();
      mat Sj = symm(Sigma.slice(j));
      Sj.diag() += eps;
      
      double dij = std::sqrt(std::max(0.0, gauss_w2_sq(mi, Si, mj, Sj, eps)));
      D(i,j) = dij;
      D(j,i) = dij;
    }
  }
  
  // Initialize U on Stiefel
  mat U = randn<mat>(p, (uword)d);
  U = qr_retraction(U);
  
  vec obj_hist(max_iter, fill::zeros);
  vec grad_hist(max_iter, fill::zeros);
  vec step_hist(max_iter, fill::zeros);
  ivec ls_hist(max_iter, fill::zeros);
  
  double obj_prev = datum::inf;
  
  for (int it = 0; it < max_iter; ++it) {
    // Compute objective and Euclidean gradient at current U
    double obj = 0.0;
    mat gradE(p, (uword)d, fill::zeros);
    
    for (uword i = 0; i < n; ++i) {
      vec mi = M.row(i).t();
      mat Si = symm(Sigma.slice(i));
      Si.diag() += eps;
      
      for (uword j = i+1; j < n; ++j) {
        vec mj = M.row(j).t();
        mat Sj = symm(Sigma.slice(j));
        Sj.diag() += eps;
        
        double s_ij;
        mat grad_s;
        pair_s_and_gradU(U, mi, Si, mj, Sj, s_ij, grad_s, eps);
        
        double delta = std::sqrt(std::max(s_ij, 0.0) + eps);
        double Dij   = D(i,j);
        double r = (delta - Dij);
        
        obj += r * r;
        
        // ∇L contribution: ((δ - D)/δ) * ∇s
        gradE += (r / delta) * grad_s;
      }
    }
    
    // Riemannian gradient
    mat UtG = U.t() * gradE;
    mat gradR = gradE - U * symm(UtG);
    double gnorm = norm(gradR, "fro");
    
    obj_hist(it) = obj;
    grad_hist(it) = gnorm;
    
    if (verbose) {
      Rcpp::Rcout << "iter " << it+1 << "  obj=" << obj << "  |grad|=" << gnorm;
    }
    
    // stopping
    if (gnorm < tol_grad) {
      if (verbose) Rcpp::Rcout << "  (stop: tol_grad)\n";
      obj_hist = obj_hist.head(it+1);
      grad_hist = grad_hist.head(it+1);
      step_hist = step_hist.head(it+1);
      ls_hist   = ls_hist.head(it+1);
      break;
    }
    if (std::abs(obj_prev - obj) < tol_obj * (1.0 + std::abs(obj_prev))) {
      if (verbose) Rcpp::Rcout << "  (stop: tol_obj)\n";
      obj_hist = obj_hist.head(it+1);
      grad_hist = grad_hist.head(it+1);
      step_hist = step_hist.head(it+1);
      ls_hist   = ls_hist.head(it+1);
      break;
    }
    obj_prev = obj;
    
    // ---- Armijo backtracking line search on Stiefel ----
    double alpha = step_init;
    double g2 = gnorm * gnorm;
    
    bool accepted = false;
    int ls_used = 0;
    
    mat U_new;
    
    for (int ls = 0; ls < ls_max; ++ls) {
      // candidate via QR retraction
      U_new = qr_retraction(U - alpha * gradR);
      
      double obj_new = objective_only(U_new, M, Sigma, D, eps);
      
      // Armijo sufficient decrease:
      // L(U_new) <= L(U) - c * alpha * ||grad||^2
      if (obj_new <= obj - armijo_c * alpha * g2) {
        accepted = true;
        ls_used = ls + 1;
        obj_prev = obj; // keep previous for tol_obj check next iter
        obj = obj_new;  // not strictly necessary, but useful if you print
        break;
      }
      
      alpha *= armijo_beta;
      if (alpha < ls_min_step) break;
    }
    
    step_hist(it) = alpha;
    ls_hist(it)   = ls_used;
    
    if (!accepted) {
      // Could not find sufficient decrease (typically numerical / very flat region).
      // You can either break or take the smallest step; breaking is safer.
      if (verbose) Rcpp::Rcout << "  (line search failed)\n";
      obj_hist = obj_hist.head(it+1);
      grad_hist = grad_hist.head(it+1);
      step_hist = step_hist.head(it+1);
      ls_hist   = ls_hist.head(it+1);
      break;
    }
    
    if (verbose) Rcpp::Rcout << "  step=" << alpha << "  ls=" << ls_used << "\n";
    
    // accept
    U = U_new;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("U") = U,
    Rcpp::Named("D") = D,
    Rcpp::Named("obj") = obj_hist,
    Rcpp::Named("grad_norm") = grad_hist,
    Rcpp::Named("step") = step_hist,
    Rcpp::Named("ls_iters") = ls_hist
  );
}
