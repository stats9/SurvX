// survcore: small numerical summaries shared by the front ends (distribution functions and their
// inverses, sample quantiles, confidence intervals on transformed scales, frailty summaries,
// Nelson-Aalen tables, step functions). Everything a front end would otherwise recompute is here,
// so that R and other front ends give the same numbers.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "special.hpp"

namespace survcore {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// ---------------------------------------------------------------------------------------------
// distributions
// ---------------------------------------------------------------------------------------------
// regularised lower incomplete gamma P(a, x), accurate in both tails (series for x < a + 1,
// otherwise 1 - Q from the continued fraction of gamma_q)
inline double gamma_p(double a, double x) {
  if (!(x > 0.0)) return 0.0;
  if (!(a > 0.0)) return kNaN;
  if (std::isinf(x)) return 1.0;
  if (x < a + 1.0) {
    double ap = a, sum = 1.0 / a, del = sum;
    for (int n = 0; n < 10000; ++n) {
      ap += 1.0; del *= x / ap; sum += del;
      if (std::fabs(del) < std::fabs(sum) * 1e-17) break;
    }
    return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
  }
  return 1.0 - gamma_q(a, x);
}

// standard normal distribution function and its inverse (Acklam's approximation refined by one
// Halley step: full double precision; Acklam 2003)
inline double pnorm_std(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
inline double qnorm_p(double p) {
  if (!(p > 0.0)) return p == 0.0 ? -INFINITY : kNaN;
  if (!(p < 1.0)) return p == 1.0 ? INFINITY : kNaN;
  if (p > 0.5) return -qnorm_p(1.0 - p);                 // 1 - p is exact for p >= 0.5 (Sterbenz)
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
  double x;
  if (p < 0.02425) { const double q = std::sqrt(-2 * std::log(p)); x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  else if (p > 1 - 0.02425) { const double q = std::sqrt(-2 * std::log1p(-p)); x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  else { const double q = p - 0.5, r = q * q; x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1); }
  const double e = pnorm_std(x) - p;
  const double u = e * 2.5066282746310002 * std::exp(x * x / 2);
  x -= u / (1 + x * u / 2);
  return x;
}

// quantile of a gamma distribution (shape a, rate r): Newton on y = log x, safeguarded by
// bisection on y; start from the small-x expansion P(a, x) ~ x^a / Gamma(a + 1) when it applies,
// otherwise from the Wilson-Hilferty approximation
inline double qgamma_rate(double p, double a, double r) {
  if (!(a > 0.0) || !(r > 0.0) || !(p >= 0.0 && p <= 1.0)) return kNaN;
  if (p == 0.0) return 0.0;
  if (p == 1.0) return INFINITY;
  const double lg = std::lgamma(a);
  double y = (std::log(p) + std::lgamma(a + 1.0)) / a;       // log x of the small-x expansion
  if (!(y < -1.0)) {
    const double z = qnorm_p(p), h = 1.0 / (9.0 * a);
    const double c = 1.0 - h + z * std::sqrt(h);
    y = c > 0.0 ? std::log(a) + 3.0 * std::log(c) : std::log(a) - 3.0;
  }
  double lo = -INFINITY, hi = INFINITY;
  for (int it = 0; it < 300; ++it) {
    const double x = std::exp(y);
    const double f = gamma_p(a, x) - p;
    if (f == 0.0) break;
    if (f < 0.0) lo = y; else hi = y;
    const double dfdy = std::exp(a * y - x - lg);             // d P / d log x = x^a e^-x / Gamma(a)
    double yn = dfdy > 0.0 ? y - f / dfdy : kNaN;
    if (!(yn > lo && yn < hi) || !std::isfinite(yn)) {
      if (std::isfinite(lo) && std::isfinite(hi)) yn = 0.5 * (lo + hi);
      else yn = std::isfinite(lo) ? y + 2.0 : y - 2.0 * std::fmax(1.0, std::fabs(y));
    }
    if (std::fabs(yn - y) <= 1e-15 * std::fmax(1.0, std::fabs(y))) { y = yn; break; }
    y = yn;
  }
  return std::exp(y) / r;
}

// upper tail of the chi-square distribution
inline double pchisq_upper(double x, double df) {
  if (std::isnan(x) || !(df > 0.0)) return kNaN;
  if (!(x > 0.0)) return 1.0;
  if (std::isinf(x)) return 0.0;
  return gamma_q(0.5 * df, 0.5 * x);
}

// ---------------------------------------------------------------------------------------------
// sample quantiles of Hyndman & Fan (1996), computed exactly as R's quantile(); x sorted and
// finite, p in [0, 1]
// ---------------------------------------------------------------------------------------------
inline double quantile7_sorted(const std::vector<double>& x, double p) {
  const std::size_t n = x.size();
  if (!n) return kNaN;
  const double index = 1.0 + static_cast<double>(n - 1) * p;   // as R's quantile(type = 7)
  const double lo = std::floor(index), hi = std::ceil(index);
  const double qlo = x[static_cast<std::size_t>(lo) - 1], qhi = x[static_cast<std::size_t>(hi) - 1];
  if (index > lo && qhi != qlo) { const double h = index - lo; return (1.0 - h) * qlo + h * qhi; }
  return qlo;
}

// type 1 (inverse of the empirical distribution function; R's quantile(type = 1))
inline double quantile1_sorted(const std::vector<double>& x, double p) {
  const std::size_t n = x.size();
  if (!n) return kNaN;
  const double np = static_cast<double>(n) * p;
  std::size_t j = static_cast<std::size_t>(std::floor(np + 4.0 * std::numeric_limits<double>::epsilon()));
  if (np > static_cast<double>(j)) ++j;                           // as R's quantile(type = 1), 1-based
  if (j < 1) j = 1;
  if (j > n) j = n;
  return x[j - 1];
}

// ---------------------------------------------------------------------------------------------
// confidence intervals
// ---------------------------------------------------------------------------------------------
// cumulative incidence F with standard error se on the complementary log-log scale:
// F^exp(-/+ z se / (F log F)); NaN where F is 0 or 1 or se is not finite
inline void cif_interval(double F, double se, double zq, double& lower, double& upper) {
  if (!(std::isfinite(F) && std::isfinite(se) && F > 0.0 && F < 1.0)) { lower = upper = kNaN; return; }
  const double e = zq * se / (F * std::log(F));
  lower = std::pow(F, std::exp(-e));
  upper = std::pow(F, std::exp(e));
}

// ---------------------------------------------------------------------------------------------
// step functions
// ---------------------------------------------------------------------------------------------
// value of a right-continuous step function with jumps at increasing times t[0..m-1] (value v[j]
// from t[j] on) at time s; `before` before the first jump
inline double step_at(const double* t, const double* v, std::size_t m, double s, double before) {
  const double* it = std::upper_bound(t, t + m, s);
  return it == t ? before : v[(it - t) - 1];
}

// ---------------------------------------------------------------------------------------------
// shared gamma frailty
// ---------------------------------------------------------------------------------------------
// likelihood ratio test of no frailty: half chi-square with 1 df (Self & Liang 1987)
inline void frailty_lrt(double ll, double ll0, bool boundary, double& stat, double& p) {
  stat = boundary ? 0.0 : std::fmax(0.0, 2.0 * (ll - ll0));
  p = 0.5 * pchisq_upper(stat, 1.0);
}

}  // namespace survcore
