// survcore: scaling of covariates for the fits, and detection of collinear covariates.
//
// Column scaling. A covariate on a very large or very small scale (income in rials, time in
// seconds) makes Newton-Raphson badly conditioned. Every fit therefore works with X / s, s a power
// of two chosen per column so that the root mean square of x / s lies in [1/16, 16]; columns already
// in that range keep s = 1. Because s is a power of two, x / s and the back-transformed coefficients
// b = b* / s and covariances V = V* / (s s') are exact: for data on a sensible scale nothing changes
// at all, and for the others the fit is that of the well-scaled problem.
//
// Collinearity. collinear_columns() finds covariates that are (numerically) linear combinations of
// the others and of the constant: Householder QR with column pivoting of the centred columns, each
// scaled to unit length, with tolerance 1e-7 (as lm()). The columns returned are the ones that come
// last in the pivot order, so that the earlier of two collinear columns is kept.
#ifndef SURVCORE_SCALING_HPP
#define SURVCORE_SCALING_HPP

#include <cmath>
#include <cstddef>
#include <vector>

namespace survcore {

struct ColumnScale {
  std::vector<double> s;          // scale of every column (powers of two; 1 for most)
  std::vector<double> X;          // the scaled matrix, row-major (empty when no column is scaled)
  bool any = false;
  const double* data(const double* X0) const { return any ? X.data() : X0; }
};

// the power of two that brings the root mean square of a column into [1/16, 16] (1 when it is there)
inline double column_scale(std::size_t n, std::size_t p, const double* X, std::size_t j) {
  long double ss = 0.0L;
  std::size_t k = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const double v = X[i * p + j];
    if (std::isfinite(v)) { ss += static_cast<long double>(v) * v; ++k; }
  }
  if (k == 0) return 1.0;
  const double rms = std::sqrt(static_cast<double>(ss / static_cast<long double>(k)));
  if (!(rms > 0) || !std::isfinite(rms) || (rms >= 1.0 / 16 && rms <= 16)) return 1.0;
  int e = 0;
  std::frexp(rms, &e);                                    // rms = f 2^e, f in [0.5, 1)
  return std::ldexp(1.0, e);
}

inline ColumnScale scale_columns(std::size_t n, std::size_t p, const double* X) {
  ColumnScale cs;
  cs.s.assign(p, 1.0);
  if (!X) return cs;
  for (std::size_t j = 0; j < p; ++j) { cs.s[j] = column_scale(n, p, X, j); if (cs.s[j] != 1.0) cs.any = true; }
  if (!cs.any) return cs;
  cs.X.assign(X, X + n * p);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < p; ++j) cs.X[i * p + j] /= cs.s[j];
  return cs;
}

// coefficients at positions pos[0..p) of a parameter vector: to the scaled problem (b s) and back (b / s)
inline void coef_to_scaled(const ColumnScale& cs, std::size_t first, double* theta) {
  if (cs.any) for (std::size_t j = 0; j < cs.s.size(); ++j) theta[first + j] *= cs.s[j];
}
inline void coef_from_scaled(const ColumnScale& cs, std::size_t first, double* theta) {
  if (cs.any) for (std::size_t j = 0; j < cs.s.size(); ++j) theta[first + j] /= cs.s[j];
}
// covariance P x P (row-major) of a parameter vector whose covariate coefficients start at first
inline void var_from_scaled(const ColumnScale& cs, std::size_t first, std::size_t P, double* V) {
  if (!cs.any || !V) return;
  std::vector<double> f(P, 1.0);
  for (std::size_t j = 0; j < cs.s.size(); ++j) f[first + j] = 1.0 / cs.s[j];
  for (std::size_t a = 0; a < P; ++a)
    for (std::size_t b = 0; b < P; ++b) V[a * P + b] *= f[a] * f[b];
}

// covariates (0-based) that are linear combinations of the others and of the constant
inline std::vector<int> collinear_columns(std::size_t n, std::size_t p, const double* X, double tol = 1e-7) {
  std::vector<int> dep;
  if (p == 0 || n == 0) return dep;
  std::vector<double> A(n * p);                           // column-major centred, unit-length columns
  std::vector<double> norm0(p, 0.0);
  for (std::size_t j = 0; j < p; ++j) {
    long double m = 0.0L;
    for (std::size_t i = 0; i < n; ++i) m += X[i * p + j];
    m /= static_cast<long double>(n);
    double ss = 0.0;
    for (std::size_t i = 0; i < n; ++i) { const double v = static_cast<double>(X[i * p + j] - m); A[j * n + i] = v; ss += v * v; }
    norm0[j] = std::sqrt(ss);
    if (norm0[j] > 0) for (std::size_t i = 0; i < n; ++i) A[j * n + i] /= norm0[j];
  }
  std::vector<std::size_t> piv(p);
  for (std::size_t j = 0; j < p; ++j) piv[j] = j;
  std::vector<double> rest(p);
  for (std::size_t j = 0; j < p; ++j) rest[j] = norm0[j] > 0 ? 1.0 : 0.0;
  std::size_t rank = 0;
  for (std::size_t k = 0; k < p && k < n; ++k) {
    // the remaining column with the largest residual norm
    std::size_t best = k;
    for (std::size_t j = k; j < p; ++j) if (rest[piv[j]] > rest[piv[best]]) best = j;
    std::swap(piv[k], piv[best]);
    const std::size_t c = piv[k];
    double nrm = 0.0;
    for (std::size_t i = k; i < n; ++i) nrm += A[c * n + i] * A[c * n + i];
    nrm = std::sqrt(nrm);
    if (!(nrm > tol)) break;                              // every remaining column is dependent
    ++rank;
    const double alpha = A[c * n + k] > 0 ? -nrm : nrm;
    std::vector<double> u(n - k);
    for (std::size_t i = k; i < n; ++i) u[i - k] = A[c * n + i];
    u[0] -= alpha;
    double uu = 0.0;
    for (double v : u) uu += v * v;
    if (uu > 0)
      for (std::size_t jj = k; jj < p; ++jj) {
        const std::size_t d = piv[jj];
        double t = 0.0;
        for (std::size_t i = k; i < n; ++i) t += u[i - k] * A[d * n + i];
        t = 2.0 * t / uu;
        for (std::size_t i = k; i < n; ++i) A[d * n + i] -= t * u[i - k];
      }
    for (std::size_t jj = k + 1; jj < p; ++jj) {
      const std::size_t d = piv[jj];
      double r = 0.0;
      for (std::size_t i = k + 1; i < n; ++i) r += A[d * n + i] * A[d * n + i];
      rest[d] = std::sqrt(r);
    }
  }
  std::vector<bool> is_dep(p, false);
  for (std::size_t k = rank; k < p; ++k) is_dep[piv[k]] = true;
  for (std::size_t j = 0; j < p; ++j) if (is_dep[j]) dep.push_back(static_cast<int>(j));
  return dep;
}

}  // namespace survcore

#endif
