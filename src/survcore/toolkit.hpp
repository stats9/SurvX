// survcore: the small computations on results that front ends would otherwise write themselves:
// elementwise transforms, reductions (sum, mean, ...), column and group means, running sums,
// linear predictors and proportional-hazards predictions, curves from event counts, information
// criteria and likelihood ratio tests, Wald tables, standard errors from a covariance matrix and
// scaled Schoenfeld residuals. Means and sums follow R (long double accumulation, the correction
// pass of mean()), so that the R front end gives the numbers of R's own functions.
#ifndef SURVCORE_TOOLKIT_HPP
#define SURVCORE_TOOLKIT_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "summaries.hpp"

namespace survcore {

enum EwCode : int {
  EW_EXP = 1, EW_LOG = 2, EW_PLOGIS = 3, EW_LOGIT = 4, EW_QNORM = 5, EW_COMPLEMENT = 6, EW_NEGLOG = 7,
  EW_SUB = 8, EW_DIV = 9, EW_LOGRATIO = 10, EW_PVALUE = 11, EW_EXP_SE = 12, EW_IC_MID = 13, EW_NEGLOG_FLOOR = 14,
  EW_MUL = 15, EW_SQRT = 16, EW_ADD = 17, EW_QNORM_COMP = 18, EW_LOGIT_COMP = 19, EW_EXP_NEG = 20, EW_ABSDIFF = 21, EW_SIGNIF = 22, EW_ROUND = 23
};
inline bool ew_binary(int code) {
  return code == EW_SUB || code == EW_DIV || code == EW_LOGRATIO || code == EW_PVALUE || code == EW_EXP_SE || code == EW_IC_MID ||
         code == EW_MUL || code == EW_ADD || code == EW_ABSDIFF || code == EW_SIGNIF;
}
inline bool ew_valid(int code) { return code >= EW_EXP && code <= EW_ROUND; }

inline double plogis1(double x) {                         // as R: no overflow on either side
  if (std::isnan(x)) return x;
  if (x >= 0) return 1.0 / (1.0 + std::exp(-x));
  const double e = std::exp(x);
  return e / (1.0 + e);
}

// R_pow_di and fprec of R (signif): x rounded to the given number of significant digits
inline double r_pow_di(double x, int n) {
  double xn = 1.0;
  if (n < 0) { n = -n; x = 1 / x; }
  for (;;) { if (n & 1) xn *= x; if (n >>= 1) x *= x; else break; }
  return xn;
}
inline double r_signif(double x, double digits) {
  if (std::isnan(x) || std::isnan(digits)) return x + digits;
  if (!std::isfinite(x)) return x;
  if (!std::isfinite(digits)) { if (digits > 0.0) return x; digits = 1.0; }
  if (x == 0) return x;
  int dig = static_cast<int>(std::round(digits));
  if (dig > 22) return x;
  if (dig < 1) dig = 1;
  double sgn = 1.0;
  if (x < 0.0) { sgn = -sgn; x = -x; }
  const double l10 = std::log10(x);
  int e10 = static_cast<int>(dig - 1 - std::floor(l10));
  const int max10e = 308;
  if (std::fabs(l10) < max10e - 2) {
    double p10 = 1.0;
    if (e10 > max10e) { p10 = r_pow_di(10., e10 - max10e); e10 = max10e; }
    if (e10 > 0) { const double pw = r_pow_di(10., e10); return sgn * (std::nearbyint((x * pw) * p10) / pw) / p10; }
    const double pw = r_pow_di(10., -e10);
    return sgn * (std::nearbyint(x / pw) * pw);
  }
  // very large or small: as R, in two steps through the decimal exponent
  const bool do_round = max10e - l10 >= r_pow_di(10., -dig);
  const int e2 = dig + ((e10 > 0) ? 1 : 6);
  const double p10 = r_pow_di(10., e2);
  const double P10 = r_pow_di(10., e10 - e2);
  x *= p10; x *= P10;
  if (do_round) x += 0.5;
  x = std::floor(x) / p10;
  return sgn * x / P10;
}

inline double ew1(int code, double x, double y) {
  switch (code) {
    case EW_EXP: return std::exp(x);
    case EW_LOG: return std::log(x);
    case EW_PLOGIS: return plogis1(x);
    case EW_LOGIT: return std::log(x / (1.0 - x));
    case EW_QNORM: return std::isnan(x) ? x : (x < 0 || x > 1 ? std::nan("") : qnorm_p(x));
    case EW_COMPLEMENT: return 1.0 - x;
    case EW_NEGLOG: return -std::log(x);
    case EW_SUB: return x - y;
    case EW_DIV: return x / y;
    case EW_LOGRATIO: return std::log(x / y);
    case EW_PVALUE: { const double r = x / y; return std::isnan(r) ? r : 2.0 * pnorm_std(-std::fabs(r)); }
    case EW_EXP_SE: return std::exp(x) * std::sqrt(y);
    case EW_IC_MID: return !std::isfinite(y) ? x : (x == 0 ? y / 2.0 : (x + y) / 2.0);
    case EW_NEGLOG_FLOOR: return std::isnan(x) ? x : -std::log(std::fmax(x, 1e-300));
    case EW_MUL: return x * y;
    case EW_SQRT: return std::sqrt(x);
    case EW_ADD: return x + y;
    case EW_QNORM_COMP: { const double c = 1.0 - x; return std::isnan(c) ? c : (c < 0 || c > 1 ? std::nan("") : qnorm_p(c)); }
    case EW_LOGIT_COMP: return std::log((1.0 - x) / x);
    case EW_ABSDIFF: return std::fabs(x - y);
    case EW_SIGNIF: return r_signif(x, y);
    case EW_ROUND: return std::nearbyint(x);              // round(x) of R: halves to even
    default: return std::exp(-x);                         // EW_EXP_NEG
  }
}

// ---- reductions as R: sum (long double), mean (long double with the correction pass) ----------
enum ReduceCode : int { RD_SUM = 0, RD_MEAN = 1, RD_MAXABS = 2, RD_MAX = 3, RD_MIN = 4 };

inline double r_sum(const std::vector<double>& x) {
  long double s = 0; for (double v : x) s += v; return static_cast<double>(s);
}
inline double r_mean_of(const std::vector<double>& x) {
  if (x.empty()) return std::nan("");
  const long double n = static_cast<long double>(x.size());
  long double s = 0; for (double v : x) s += v;
  s /= n;
  if (!std::isfinite(static_cast<double>(s))) return static_cast<double>(s);
  long double t = 0; for (double v : x) t += v - s;
  return static_cast<double>(s + t / n);
}
// na_rm: leave out NaN (R's NA and NaN); otherwise a NaN gives NaN
inline double reduce(int code, std::size_t n, const double* x, std::size_t stride, bool na_rm) {
  std::vector<double> v;
  for (std::size_t i = 0; i < n; ++i) {
    const double y = x[i * stride];
    if (std::isnan(y)) { if (na_rm) continue; return y; }
    v.push_back(y);
  }
  switch (code) {
    case RD_SUM: return r_sum(v);
    case RD_MEAN: return r_mean_of(v);
    default: {
      if (v.empty()) return code == RD_MIN ? INFINITY : -INFINITY;     // as R's max/min of nothing
      double m = code == RD_MAXABS ? std::fabs(v[0]) : v[0];
      for (double y : v) {
        const double a = code == RD_MAXABS ? std::fabs(y) : y;
        if (code == RD_MIN ? a < m : a > m) m = a;
      }
      return m;
    }
  }
}

// ---- curves from event counts ---------------------------------------------------------------------
// mode 0 Nelson-Aalen H = cumsum(d / Y), 1 exp(-H), 2 Kaplan-Meier cumprod(1 - d / Y); running sums
// and products in long double, as R's cumsum() and cumprod()
inline void np_curve_inc(std::size_t m, const double* d, const double* Y, int mode, double* out) {
  long double H = 0.0L, S = 1.0L;
  for (std::size_t j = 0; j < m; ++j) {
    const double h = d[j] / Y[j];
    H += h; S *= static_cast<long double>(1.0 - h);
    out[j] = mode == 0 ? static_cast<double>(H) : (mode == 1 ? std::exp(-static_cast<double>(H)) : static_cast<double>(S));
  }
}

// running sum along each row of M (nr x nc, row-major) as R's cumsum(), then capped at cap (cap < 0: none)
inline void cumsum_rows(std::size_t nr, std::size_t nc, const double* M, double cap, double* out) {
  for (std::size_t i = 0; i < nr; ++i) {
    long double s = 0.0L;
    for (std::size_t j = 0; j < nc; ++j) {
      s += M[i * nc + j];
      const double v = static_cast<double>(s);
      out[i * nc + j] = cap >= 0 && v > cap ? cap : v;
    }
  }
}


// ---- grids ----------------------------------------------------------------------------------------
// seq(from, to, length.out = n) as R's seq.default: c(from, from + (1:(n - 2)) * ((to - from) / (n - 1)), to)
inline std::vector<double> seq_length(double from, double to, std::size_t n) {
  std::vector<double> out;
  if (n == 0) return out;
  if (n == 1) { out.push_back(from); return out; }
  if (n == 2) { out.push_back(from); out.push_back(to); return out; }
  if (from == to) return std::vector<double>(n, from);
  const double by = (to - from) / static_cast<double>(n - 1);
  out.push_back(from);
  for (std::size_t i = 1; i + 1 < n; ++i) out.push_back(from + static_cast<double>(i) * by);
  out.push_back(to);
  return out;
}

// seq(from, to, by) as R: from + (0:k) * by, k = floor((to - from) / by + 1e-10), kept within to
inline std::vector<double> seq_by(double from, double to, double by) {
  std::vector<double> out;
  const double del = to - from;
  if (del == 0) { out.push_back(from); return out; }
  const double n = del / by;
  if (!(n >= 0) || !std::isfinite(n)) return out;
  const std::size_t k = static_cast<std::size_t>(n + 1e-10);
  for (std::size_t i = 0; i <= k; ++i) {
    double v = from + static_cast<double>(i) * by;
    if (by > 0 ? v > to : v < to) v = to;
    out.push_back(v);
  }
  return out;
}

// pretty(c(lo, up), n) as R's pretty.default (R_pretty with min.n = n %/% 3, shrink.sml = 0.75,
// high.u.bias = 1.5, u5.bias = 0.5 + 1.5 * 1.5, eps.correct = 0, bounds = TRUE)
inline std::vector<double> pretty_r(double lo, double up, int ndiv) {
  const double h = 1.5, h5 = 0.5 + 1.5 * h, f_min = std::ldexp(1.0, -20), rounding_eps = 1e-10, shrink = 0.75;
  const int min_n = ndiv / 3;
  const double dx = up - lo;
  double cell;
  bool i_small;
  if (dx == 0 && up == 0) { cell = 1; i_small = true; }
  else {
    cell = std::fmax(std::fabs(lo), std::fabs(up));
    double U = 1 + ((h5 >= 1.5 * h + .5) ? 1 / (1 + h) : 1.5 / (1 + h5));
    U *= std::max(1, ndiv) * 2.220446049250313e-16;
    i_small = dx < cell * U * 3;
  }
  if (i_small) {
    if (cell > 10) cell = 9 + cell / 10;
    cell *= shrink;
    if (min_n > 1) cell /= min_n;
  } else {
    cell = dx;
    if (ndiv > 1) cell /= ndiv;
  }
  const double subsmall = f_min * 2.2250738585072014e-308;
  if (cell < subsmall) cell = subsmall;
  const double base = std::pow(10.0, std::floor(std::log10(cell)));
  double unit = base, ns;
  if ((ns = 2 * base) - cell < h * (cell - unit)) {
    unit = ns;
    if ((ns = 5 * base) - cell < h5 * (cell - unit)) {
      unit = ns;
      if ((ns = 10 * base) - cell < h * (cell - unit)) unit = ns;
    }
  }
  ns = std::floor(lo / unit + rounding_eps);
  double nu = std::ceil(up / unit - rounding_eps);
  while (ns * unit > lo + rounding_eps * unit) ns--;
  while (nu * unit < up - rounding_eps * unit) nu++;
  int k = static_cast<int>(0.5 + nu - ns);
  int nd;
  if (k < min_n) {
    k = min_n - k;
    if (ns >= 0.) { nu += k / 2; ns -= k / 2 + k % 2; }
    else { ns -= k / 2; nu += k / 2 + k % 2; }
    nd = min_n;
  } else nd = k;
  std::vector<double> s = seq_by(ns * unit, nu * unit, unit);
  if (nd > 0) {
    const double delta = (nu * unit - ns * unit) / nd;
    for (double& v : s) if (std::fabs(v) < 1e-14 * delta) v = 0.0;
  }
  return s;
}

}  // namespace survcore

#endif
