// survcore: minimal dense linear algebra for small (p x p) symmetric systems.
// Deliberately dependency-free. Reference for the Cholesky factorisation:
//   Golub & Van Loan (2013), Matrix Computations, 4th ed., Algorithm 4.2.1.
#pragma once
#include <cmath>
#include <cstddef>
#include <vector>

namespace survcore {

struct Mat {
  std::size_t n = 0;            // square n x n, row-major
  std::vector<double> a;
  Mat() = default;
  explicit Mat(std::size_t n_) : n(n_), a(n_ * n_, 0.0) {}
  double& operator()(std::size_t i, std::size_t j) { return a[i * n + j]; }
  double operator()(std::size_t i, std::size_t j) const { return a[i * n + j]; }
  void zero() { std::fill(a.begin(), a.end(), 0.0); }
};

// In-place lower Cholesky A = L L^T. Returns false when A is not numerically
// positive definite (pivot <= rel_tol * largest diagonal). The caller must
// report that honestly instead of "repairing" the matrix.
inline bool cholesky(Mat& A, double rel_tol = 1e-12) {
  const std::size_t n = A.n;
  double dmax = 0.0;
  for (std::size_t i = 0; i < n; ++i) dmax = std::fmax(dmax, std::fabs(A(i, i)));
  if (!(dmax > 0.0) || !std::isfinite(dmax)) return false;
  for (std::size_t j = 0; j < n; ++j) {
    double d = A(j, j);
    for (std::size_t k = 0; k < j; ++k) d -= A(j, k) * A(j, k);
    if (!(d > rel_tol * dmax) || !std::isfinite(d)) return false;
    d = std::sqrt(d);
    A(j, j) = d;
    for (std::size_t i = j + 1; i < n; ++i) {
      double s = A(i, j);
      for (std::size_t k = 0; k < j; ++k) s -= A(i, k) * A(j, k);
      A(i, j) = s / d;
    }
  }
  return true;
}

// Solve L L^T x = b given the factor from cholesky().
inline void chol_solve(const Mat& L, std::vector<double>& b) {
  const std::size_t n = L.n;
  for (std::size_t i = 0; i < n; ++i) {
    double s = b[i];
    for (std::size_t k = 0; k < i; ++k) s -= L(i, k) * b[k];
    b[i] = s / L(i, i);
  }
  for (std::size_t ii = n; ii-- > 0;) {
    double s = b[ii];
    for (std::size_t k = ii + 1; k < n; ++k) s -= L(k, ii) * b[k];
    b[ii] = s / L(ii, ii);
  }
}

// Inverse of the original matrix from its Cholesky factor.
inline Mat chol_inverse(const Mat& L) {
  const std::size_t n = L.n;
  Mat inv(n);
  std::vector<double> e(n);
  for (std::size_t j = 0; j < n; ++j) {
    std::fill(e.begin(), e.end(), 0.0);
    e[j] = 1.0;
    chol_solve(L, e);
    for (std::size_t i = 0; i < n; ++i) inv(i, j) = e[i];
  }
  return inv;
}

}  // namespace survcore
