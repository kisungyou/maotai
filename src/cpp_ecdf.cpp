// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// ---------- helpers -------------------------------------------------------

// sup_k |a_k - b_k|
static inline double max_abs_diff(const vec& a, const vec& b) {
  double m = 0.0;
  const uword L = a.n_elem;
  for (uword k = 0; k < L; ++k) {
    double v = std::fabs(a[k] - b[k]);
    if (v > m) m = v;
  }
  return m;
}

// Trapezoid integral of |a-b|^p over irregular grid with increments d (length L-1).
// If use_inf, return sup norm instead.
static inline double lp_trap(const vec& a,
                             const vec& b,
                             double p,
                             const vec& d,
                             bool use_inf) {
  const uword L = a.n_elem;
  if (use_inf) return max_abs_diff(a, b);
  
  double acc = 0.0;
  double g_prev = std::pow(std::fabs(a[0] - b[0]), p);
  for (uword k = 0; k + 1 < L; ++k) {
    double g_curr = std::pow(std::fabs(a[k+1] - b[k+1]), p);
    acc += d[k] * (g_prev + g_curr) * 0.5;
    g_prev = g_curr;
  }
  return acc; // caller raises to 1/p
}

// ---------- exported kernels ----------------------------------------------

// F: (L x N) columns are ECDF values on a shared grid (from your elist_fform()).
// Returns N x N symmetric matrix of KS distances.
// [[Rcpp::export]]
arma::mat cpp_dist_ks(const arma::mat& F) {
  const uword L = F.n_rows;
  const uword N = F.n_cols;
  arma::mat D(N, N, fill::zeros);
  
  vec fi(L), fj(L);
  for (uword i = 0; i < N; ++i) {
    fi = F.col(i);
    for (uword j = i + 1; j < N; ++j) {
      fj = F.col(j);
      double dij = max_abs_diff(fi, fj);
      D(i,j) = dij;
      D(j,i) = dij;
    }
  }
  return D;
}

// t: grid (length L), F: (L x N) ECDF values on t, p in (0,∞].
// Returns N x N symmetric matrix of Lp distances between ECDFs over t.
// [[Rcpp::export]]
arma::mat cpp_dist_lp(const arma::vec& t, const arma::mat& F, double p) {
  const uword L = F.n_rows;
  const uword N = F.n_cols;
  arma::mat D(N, N, fill::zeros);
  
  // precompute dt
  arma::vec dt( (L > 0) ? (L - 1) : 0 );
  for (uword k = 0; k + 1 < L; ++k) dt[k] = t[k+1] - t[k];
  
  const bool use_inf = !R_finite(p);
  
  vec fi(L), fj(L);
  for (uword i = 0; i < N; ++i) {
    fi = F.col(i);
    for (uword j = i + 1; j < N; ++j) {
      fj = F.col(j);
      double val = lp_trap(fi, fj, p, dt, use_inf);
      if (!use_inf) val = std::pow(val, 1.0 / p);
      D(i,j) = val;
      D(j,i) = val;
    }
  }
  return D;
}

// q: quantile grid (length Lq), Q: (Lq x N) quantile columns (from quantile(ecdf, q)).
// Returns N x N symmetric matrix of 1D p-Wasserstein distances.
// [[Rcpp::export]]
arma::mat cpp_dist_wasserstein(const arma::vec& q,
                               const arma::mat& Q,
                               double p) {
  const uword L = Q.n_rows;
  const uword N = Q.n_cols;
  arma::mat D(N, N, fill::zeros);
  
  // precompute dq
  arma::vec dq( (L > 0) ? (L - 1) : 0 );
  for (uword k = 0; k + 1 < L; ++k) dq[k] = q[k+1] - q[k];
  
  const bool use_inf = !R_finite(p);
  
  vec Qi(L), Qj(L);
  for (uword i = 0; i < N; ++i) {
    Qi = Q.col(i);
    for (uword j = i + 1; j < N; ++j) {
      Qj = Q.col(j);
      if (use_inf) {
        double m = max_abs_diff(Qi, Qj);
        D(i,j) = m;
        D(j,i) = m;
      } else {
        double acc = 0.0;
        double g_prev = std::pow(std::fabs(Qi[0] - Qj[0]), p);
        for (uword k = 0; k + 1 < L; ++k) {
          double g_curr = std::pow(std::fabs(Qi[k+1] - Qj[k+1]), p);
          acc += dq[k] * (g_prev + g_curr) * 0.5;
          g_prev = g_curr;
        }
        double val = std::pow(acc, 1.0 / p);
        D(i,j) = val;
        D(j,i) = val;
      }
    }
  }
  return D;
}
