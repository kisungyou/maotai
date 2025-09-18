// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <algorithm>

// --- utilities ---------------------------------------------------------------

static inline arma::vec finite_sorted_copy(const arma::vec& v) {
  arma::uvec idx = arma::find_finite(v);
  arma::vec out = v.elem(idx);
  out = arma::sort(out); // ascending
  return out;
}

// Advance i over all entries equal to value 'val'; return how many advanced
static inline std::size_t advance_equal(const arma::vec& x, std::size_t& i, const double val) {
  std::size_t start = i;
  const std::size_t n = x.n_elem;
  while (i < n && x[i] == val) ++i;
  return i - start;
}

// --- core: Lp distance between ECDFs for finite p ---------------------------
// ∫ |F_x(t) - F_y(t)|^p dt  over R, using piecewise-constant difference between
// successive distinct knots in the union of {x} ∪ {y}.
//
// Important detail: on the open interval (k, next_k], the ECDF value equals
// the count of observations ≤ k divided by n (right-continuous ECDF).
//
// Complexity: sort O(n_x log n_x + n_y log n_y) + one linear merge O(n_x + n_y)
// [[Rcpp::export]]
double lp_ecdf_distance(const arma::vec& x, const arma::vec& y, double p = 2.0) {
  if (!(p > 0 && std::isfinite(p))) Rcpp::stop("p must be a positive finite number");
  arma::vec xs = finite_sorted_copy(x);
  arma::vec ys = finite_sorted_copy(y);
  const std::size_t nx = xs.n_elem, ny = ys.n_elem;
  if (nx == 0 && ny == 0) return 0.0;
  if (nx == 0 || ny == 0) Rcpp::stop("Both samples must contain at least one finite value.");
  
  std::size_t ix = 0, iy = 0;
  std::size_t cumx = 0, cumy = 0;
  double integral = 0.0;
  
  // First knot is the minimum observed value in the union
  double next_knot = std::min(xs[0], ys[0]);
  double prev_knot = next_knot;
  
  while (ix < nx || iy < ny) {
    // Current knot is min of next values in xs/ys
    double kx = (ix < nx) ? xs[ix] : std::numeric_limits<double>::infinity();
    double ky = (iy < ny) ? ys[iy] : std::numeric_limits<double>::infinity();
    double knot = std::min(kx, ky);
    
    // At 'knot', both ECDFs jump by the number of equals at that knot
    if (kx == knot) cumx += advance_equal(xs, ix, knot);
    if (ky == knot) cumy += advance_equal(ys, iy, knot);
    
    // After processing jumps at 'knot', the ECDFs are constant on (knot, next_knot]
    // Determine the next knot (next distinct value in either sample)
    double kx_next = (ix < nx) ? xs[ix] : std::numeric_limits<double>::infinity();
    double ky_next = (iy < ny) ? ys[iy] : std::numeric_limits<double>::infinity();
    next_knot = std::min(kx_next, ky_next);
    
    // Width of the interval (knot, next_knot]; if finite, accumulate contribution
    if (std::isfinite(next_knot)) {
      double Fx = static_cast<double>(cumx) / static_cast<double>(nx);
      double Fy = static_cast<double>(cumy) / static_cast<double>(ny);
      double diff = std::fabs(Fx - Fy);
      double width = next_knot - knot;
      if (width < 0) Rcpp::stop("Internal error: negative interval width.");
      if (width > 0) {
        integral += std::pow(diff, p) * width;
      }
      prev_knot = next_knot;
    } else {
      // Past the last knot, both ECDFs are 1; tail contributes zero width
      break;
    }
  }
  
  return std::pow(integral, 1.0 / p);
}

// --- KS distance (p = infinity) ---------------------------------------------
// sup_t |F_x(t) - F_y(t)|
// [[Rcpp::export]]
double ks_ecdf_distance(const arma::vec& x, const arma::vec& y) {
  arma::vec xs = finite_sorted_copy(x);
  arma::vec ys = finite_sorted_copy(y);
  const std::size_t nx = xs.n_elem, ny = ys.n_elem;
  if (nx == 0 && ny == 0) return 0.0;
  if (nx == 0 || ny == 0) Rcpp::stop("Both samples must contain at least one finite value.");
  
  std::size_t ix = 0, iy = 0;
  std::size_t cumx = 0, cumy = 0;
  double ks = 0.0;
  
  while (ix < nx || iy < ny) {
    double tx = (ix < nx) ? xs[ix] : std::numeric_limits<double>::infinity();
    double ty = (iy < ny) ? ys[iy] : std::numeric_limits<double>::infinity();
    
    if (tx <= ty) {
      double val = tx;
      while (ix < nx && xs[ix] == val) ++ix, ++cumx;
      double Fx = static_cast<double>(cumx) / nx;
      double Gy = static_cast<double>(cumy) / ny; // y not advanced yet if ty > tx
      ks = std::max(ks, std::fabs(Fx - Gy));
      if (val == ty) {
        while (iy < ny && ys[iy] == val) ++iy, ++cumy;
        Gy = static_cast<double>(cumy) / ny;
        ks = std::max(ks, std::fabs(Fx - Gy));
      }
    } else {
      double val = ty;
      while (iy < ny && ys[iy] == val) ++iy, ++cumy;
      double Fx = static_cast<double>(cumx) / nx;
      double Gy = static_cast<double>(cumy) / ny;
      ks = std::max(ks, std::fabs(Fx - Gy));
    }
  }
  return ks;
}

// --- batched: pairwise Lp distances across many samples ---------------------
// Sort each sample once; then compute pairwise distances with linear merges.
// Optionally, you can set p = INFINITY from R to request KS; otherwise finite p.
//
// [[Rcpp::export]]
arma::mat lp_ecdf_distance_matrix(const arma::field<arma::vec>& samples, double p = 2.0) {
  const std::size_t m = samples.n_elem;
  if (m == 0) return arma::mat(0,0);
  
  // Pre-sort and store finite-only copies
  std::vector<arma::vec> sorted(m);
  for (std::size_t i = 0; i < m; ++i) {
    sorted[i] = finite_sorted_copy(samples(i));
    if (sorted[i].n_elem == 0) Rcpp::stop("Sample %zu has no finite values.", i+1);
  }
  
  arma::mat D(m, m, arma::fill::zeros);
  
  auto dist_pair_finite_p = [&](const arma::vec& a, const arma::vec& b) {
    // merge-based integral using already-sorted vectors
    const std::size_t na = a.n_elem, nb = b.n_elem;
    std::size_t ia = 0, ib = 0;
    std::size_t cuma = 0, cumb = 0;
    double integral = 0.0;
    
    // first knot
    double next_knot = std::min(a[0], b[0]);
    
    while (ia < na || ib < nb) {
      double ka = (ia < na) ? a[ia] : std::numeric_limits<double>::infinity();
      double kb = (ib < nb) ? b[ib] : std::numeric_limits<double>::infinity();
      double knot = std::min(ka, kb);
      
      if (ka == knot) {
        double val = ka;
        while (ia < na && a[ia] == val) { ++ia; ++cuma; }
      }
      if (kb == knot) {
        double val = kb;
        while (ib < nb && b[ib] == val) { ++ib; ++cumb; }
      }
      
      double ka_next = (ia < na) ? a[ia] : std::numeric_limits<double>::infinity();
      double kb_next = (ib < nb) ? b[ib] : std::numeric_limits<double>::infinity();
      next_knot = std::min(ka_next, kb_next);
      
      if (std::isfinite(next_knot)) {
        double Fa = static_cast<double>(cuma) / na;
        double Fb = static_cast<double>(cumb) / nb;
        double width = next_knot - knot;
        if (width > 0) integral += std::pow(std::fabs(Fa - Fb), p) * width;
      } else {
        break;
      }
    }
    return std::pow(integral, 1.0 / p);
  };
  
  auto ks_pair = [&](const arma::vec& a, const arma::vec& b) {
    const std::size_t na = a.n_elem, nb = b.n_elem;
    std::size_t ia = 0, ib = 0;
    std::size_t cuma = 0, cumb = 0;
    double ks = 0.0;
    while (ia < na || ib < nb) {
      double ta = (ia < na) ? a[ia] : std::numeric_limits<double>::infinity();
      double tb = (ib < nb) ? b[ib] : std::numeric_limits<double>::infinity();
      if (ta <= tb) {
        double val = ta;
        while (ia < na && a[ia] == val) { ++ia; ++cuma; }
        double Fa = static_cast<double>(cuma) / na;
        double Gb = static_cast<double>(cumb) / nb;
        ks = std::max(ks, std::fabs(Fa - Gb));
        if (val == tb) {
          while (ib < nb && b[ib] == val) { ++ib; ++cumb; }
          Gb = static_cast<double>(cumb) / nb;
          ks = std::max(ks, std::fabs(Fa - Gb));
        }
      } else {
        double val = tb;
        while (ib < nb && b[ib] == val) { ++ib; ++cumb; }
        double Fa = static_cast<double>(cuma) / na;
        double Gb = static_cast<double>(cumb) / nb;
        ks = std::max(ks, std::fabs(Fa - Gb));
      }
    }
    return ks;
  };
  
  // Upper triangle
  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = i+1; j < m; ++j) {
      double dij = std::isinf(p) ? ks_pair(sorted[i], sorted[j])
        : dist_pair_finite_p(sorted[i], sorted[j]);
      D(i,j) = D(j,i) = dij;
    }
  }
  return D;
}

// --- convenience: distance to a reference sample across many -----------------
// [[Rcpp::export]]
arma::vec lp_ecdf_distance_to_reference(const arma::field<arma::vec>& samples,
                                        const arma::vec& reference,
                                        double p = 2.0) {
  const std::size_t m = samples.n_elem;
  arma::vec ref_sorted = finite_sorted_copy(reference);
  if (ref_sorted.n_elem == 0) Rcpp::stop("Reference has no finite values.");
  arma::vec out(m, arma::fill::zeros);
  
  for (std::size_t i = 0; i < m; ++i) {
    arma::vec s = finite_sorted_copy(samples(i));
    if (s.n_elem == 0) Rcpp::stop("Sample %zu has no finite values.", i+1);
    if (std::isinf(p)) {
      // KS
      out[i] = ks_ecdf_distance(s, ref_sorted);
    } else {
      out[i] = lp_ecdf_distance(s, ref_sorted, p);
    }
  }
  return out;
}
