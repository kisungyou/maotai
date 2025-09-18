// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <limits>
#include <algorithm>
#include <cmath>

// --------- small utils -------------------------------------------------------

static inline arma::uvec finite_idx(const arma::vec& v) { return arma::find_finite(v); }

static inline arma::vec finite_sorted_copy(const arma::vec& v) {
  arma::vec out = v.elem(finite_idx(v));
  if (out.n_elem == 0) Rcpp::stop("Sample has no finite values.");
  out = arma::sort(out);
  return out;
}

static inline void check_p(double p) {
  if (!(std::isfinite(p) && p >= 1.0)) Rcpp::stop("p must be finite and >= 1.");
}

static inline void normalize_weights(arma::vec& w) {
  if (w.n_elem == 0) Rcpp::stop("Empty weight vector.");
  if (arma::any(w < 0)) Rcpp::stop("Weights must be nonnegative.");
  double s = arma::accu(w);
  if (!(s > 0) || !std::isfinite(s)) Rcpp::stop("Weights must sum to a positive finite value.");
  w /= s;
}

// --------- core: Wp for unweighted samples (uniform masses) ------------------
// Returns W_p (not ^p). If sizes match, uses the order-statistics shortcut.
// Otherwise, uses a two-pointer mass-matching algorithm that transports
// masses 1/n and 1/m between sorted atoms.
// [[Rcpp::export]]
double wasserstein_p_unweighted(const arma::vec& x, const arma::vec& y, double p = 2.0) {
  check_p(p);
  arma::vec xs = finite_sorted_copy(x);
  arma::vec ys = finite_sorted_copy(y);
  const std::size_t nx = xs.n_elem, ny = ys.n_elem;
  
  if (nx == ny) {
    // Shortcut: W_p^p = (1/n) sum |x_(i) - y_(i)|^p
    arma::vec diff = arma::abs(xs - ys);
    double wp_p = arma::mean(arma::pow(diff, p));
    return std::pow(wp_p, 1.0 / p);
  }
  
  // General unequal-size case via monotone transport (two-pointer)
  std::size_t i = 0, j = 0;
  double wi = 1.0 / static_cast<double>(nx); // remaining mass at xs[i]
  double wj = 1.0 / static_cast<double>(ny); // remaining mass at ys[j]
  double wi_left = wi;
  double wj_left = wj;
  double cost_p = 0.0;
  
  while (i < nx && j < ny) {
    double move = std::min(wi_left, wj_left); // amount of mass to couple
    double gap = std::fabs(xs[i] - ys[j]);
    cost_p += move * std::pow(gap, p);
    
    wi_left -= move;
    wj_left -= move;
    
    if (wi_left == 0.0) {
      ++i;
      if (i < nx) { wi_left = wi; }
    }
    if (wj_left == 0.0) {
      ++j;
      if (j < ny) { wj_left = wj; }
    }
  }
  // By construction, both masses sum to 1, so both loops end together.
  
  return std::pow(cost_p, 1.0 / p);
}

// --------- core: Wp for weighted samples -------------------------------------
// x with weights wx, y with weights wy; weights nonnegative and sum to 1 each.
// Works for any finite p >= 1.
// [[Rcpp::export]]
double wasserstein_p_weighted(const arma::vec& x, arma::vec wx,
                              const arma::vec& y, arma::vec wy,
                              double p = 2.0) {
  check_p(p);
  if (x.n_elem == 0 || y.n_elem == 0) Rcpp::stop("Empty sample.");
  if (wx.n_elem != x.n_elem || wy.n_elem != y.n_elem) Rcpp::stop("Weight length mismatch.");
  
  normalize_weights(wx);
  normalize_weights(wy);
  
  // sort by support, permute weights accordingly, and coalesce equal locations
  auto sort_and_coalesce = [](const arma::vec& s, const arma::vec& w) {
    arma::uvec idx = finite_idx(s); // drop non-finite x; drop corresponding weights
    arma::vec sx = s.elem(idx);
    arma::vec sw = w.elem(idx);
    if (sx.n_elem == 0) Rcpp::stop("All values non-finite in one sample.");
    arma::uvec ord = arma::sort_index(sx);
    arma::vec xs = sx.elem(ord);
    arma::vec ws = sw.elem(ord);
    
    // coalesce duplicates (same coordinate) to avoid tiny steps
    std::vector<double> xv, wv;
    xv.reserve(xs.n_elem); wv.reserve(ws.n_elem);
    std::size_t i = 0, n = xs.n_elem;
    while (i < n) {
      double val = xs[i];
      double sumw = 0.0;
      while (i < n && xs[i] == val) { sumw += ws[i]; ++i; }
      xv.push_back(val);
      wv.push_back(sumw);
    }
    arma::vec xout(xv.size()), wout(wv.size());
    for (std::size_t k = 0; k < xv.size(); ++k) { xout[k] = xv[k]; wout[k] = wv[k]; }
    return std::pair<arma::vec, arma::vec>(xout, wout);
  };
  
  auto X = sort_and_coalesce(x, wx);
  auto Y = sort_and_coalesce(y, wy);
  const arma::vec& xs = X.first;   const arma::vec& ws = X.second;
  const arma::vec& ys = Y.first;   const arma::vec& wy2 = Y.second;
  
  std::size_t i = 0, j = 0;
  double wi_left = ws[0];
  double wj_left = wy2[0];
  double cost_p = 0.0;
  
  while (i < xs.n_elem && j < ys.n_elem) {
    double move = std::min(wi_left, wj_left);
    double gap = std::fabs(xs[i] - ys[j]);
    cost_p += move * std::pow(gap, p);
    
    wi_left -= move;
    wj_left -= move;
    if (wi_left <= 0) {
      ++i;
      if (i < xs.n_elem) wi_left = ws[i];
    }
    if (wj_left <= 0) {
      ++j;
      if (j < ys.n_elem) wj_left = wy2[j];
    }
  }
  // both cumulative weights are 1, so we should end together.
  return std::pow(cost_p, 1.0 / p);
}

// --------- batched pairwise distance matrix (unweighted) ---------------------
// Builds an MxM matrix of W_p distances for a list/field of numeric vectors.
// Each sample is sorted once; each pair is computed via the monotone transport.
// [[Rcpp::export]]
arma::mat wasserstein_p_distance_matrix(const arma::field<arma::vec>& samples,
                                        double p = 2.0) {
  check_p(p);
  const std::size_t m = samples.n_elem;
  if (m == 0) return arma::mat(0,0);
  
  // presort
  std::vector<arma::vec> sorted(m);
  for (std::size_t i = 0; i < m; ++i) {
    sorted[i] = finite_sorted_copy(samples(i));
  }
  
  arma::mat D(m, m, arma::fill::zeros);
  
  auto wpp_pair = [&](const arma::vec& a, const arma::vec& b) -> double {
    const std::size_t na = a.n_elem, nb = b.n_elem;
    
    if (na == nb) {
      // shortcut: (1/n sum |a-b|^p)^(1/p)
      arma::vec diff = arma::abs(a - b);
      double wp_p = arma::mean(arma::pow(diff, p));
      return std::pow(wp_p, 1.0 / p);
    }
    
    // unequal sizes: two-pointer mass matching with masses 1/na and 1/nb
    std::size_t i = 0, j = 0;
    double wi = 1.0 / static_cast<double>(na);
    double wj = 1.0 / static_cast<double>(nb);
    double wi_left = wi, wj_left = wj;
    double cost_p = 0.0;
    
    while (i < na && j < nb) {
      double move = std::min(wi_left, wj_left);
      double gap = std::fabs(a[i] - b[j]);
      cost_p += move * std::pow(gap, p);
      wi_left -= move; if (wi_left == 0.0) { ++i; if (i < na) wi_left = wi; }
      wj_left -= move; if (wj_left == 0.0) { ++j; if (j < nb) wj_left = wj; }
    }
    return std::pow(cost_p, 1.0 / p);
  };
  
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = i+1; j < m; ++j) {
      double d = wpp_pair(sorted[i], sorted[j]);
      D(i,j) = D(j,i) = d;
    }
  }
  return D;
}
