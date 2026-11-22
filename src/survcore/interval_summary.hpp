// survcore: summaries of a Turnbull NPMLE (turnbull() in R): the survival function with the mass
// of every innermost interval spread uniformly over it and the identified bounds, the curve as a
// polyline, quantiles, the weighted logrank test with Finkelstein scores and permutational
// variance, and the bootstrap of the curve.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "adjust.hpp"
#include "interval.hpp"
#include "linalg.hpp"
#include "summaries.hpp"

namespace survcore {

// P(T > t) (left = false) or P(T >= t) (left = true) from the innermost intervals (lo, hi] with
// masses p: estimate (NaN when mass of an interval reaching infinity is ambiguous) and the lower
// and upper limits over all distributions of the masses within their intervals
inline void tb_surv(std::size_t J, const double* lo, const double* hi, const double* p, double t, bool left, double& est,
                    double& lower, double& upper) {
  double sure = 0.0, amb = 0.0, spread = 0.0;
  bool undefined = false;
  for (std::size_t j = 0; j < J; ++j) {
    const bool pt = lo[j] == hi[j];
    const bool s = left ? lo[j] >= t : (pt ? lo[j] > t : lo[j] >= t);
    if (s) sure += p[j];
  }
  for (std::size_t j = 0; j < J; ++j) {
    const bool pt = lo[j] == hi[j];
    if (pt || !(lo[j] < t && hi[j] > t)) continue;
    amb += p[j];
    if (!std::isfinite(hi[j]) && p[j] > 0.0) undefined = true;
  }
  lower = sure; upper = sure + amb;
  if (undefined) { est = kNaN; return; }
  for (std::size_t j = 0; j < J; ++j) {
    const bool pt = lo[j] == hi[j];
    if (pt || !(lo[j] < t && hi[j] > t)) continue;
    spread += p[j] * (hi[j] - t) / (hi[j] - lo[j]);
  }
  est = sure + spread;
}

// the curve as a polyline: at every end point (0, the lo and the finite hi) the left limit and the
// value; non-finite points dropped
inline void tb_steps(std::size_t J, const double* lo, const double* hi, const double* p, std::vector<double>& time,
                     std::vector<double>& surv) {
  std::vector<double> g{0.0};
  for (std::size_t j = 0; j < J; ++j) { g.push_back(lo[j]); if (std::isfinite(hi[j])) g.push_back(hi[j]); }
  std::sort(g.begin(), g.end());
  g.erase(std::unique(g.begin(), g.end()), g.end());
  time.clear(); surv.clear();
  for (double t : g) {
    double e1, a, b, e2;
    tb_surv(J, lo, hi, p, t, true, e1, a, b);
    tb_surv(J, lo, hi, p, t, false, e2, a, b);
    if (std::isfinite(e1)) { time.push_back(t); surv.push_back(e1); }
    if (std::isfinite(e2)) { time.push_back(t); surv.push_back(e2); }
  }
}

// first time at which a curve on (time, surv) drops to p, linear interpolation; NaN if never
inline double tb_quantile(std::size_t m, const double* time, const double* surv, double p) {
  std::vector<double> t, s;
  for (std::size_t j = 0; j < m; ++j) if (std::isfinite(surv[j])) { t.push_back(time[j]); s.push_back(surv[j]); }
  for (std::size_t i = 0; i < s.size(); ++i) {
    if (!(s[i] <= p)) continue;
    if (i == 0) return t[0];
    if (s[i - 1] == s[i]) return t[i];
    return t[i - 1] + (s[i - 1] - p) / (s[i - 1] - s[i]) * (t[i] - t[i - 1]);
  }
  return kNaN;
}

// linear interpolation as R's approx(x, y, xout, rule = 1) (x increasing): NaN outside [x_0, x_last]
inline double approx_linear(const std::vector<double>& x, const std::vector<double>& y, double v) {
  const std::size_t n = x.size();
  if (n == 0 || std::isnan(v) || v < x[0] || v > x[n - 1]) return kNaN;
  std::size_t i = 0, j = n - 1;
  while (i < j - 1) { const std::size_t ij = (i + j) / 2; if (v < x[ij]) j = ij; else i = ij; }
  if (v == x[j]) return y[j];
  if (v == x[i]) return y[i];
  return y[i] + (y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i]));
}

// weighted logrank test for interval-censored data: Finkelstein (1986) scores from the NPMLE of the
// pooled data, permutation central limit variance (Fay & Shaw 2010). g: group codes 0..G-1.
struct TurnbullTest { double stat = kNaN, p = kNaN; int df = 0; std::vector<double> n, U, sd; };

inline TurnbullTest tb_test(int n, const double* L, const double* R, const int* g, int G, int max_iter, double tol) {
  TurnbullTest out;
  const TurnbullSupport s = turnbull_support(n, L, R);
  const TurnbullResult r = turnbull_fit(s, nullptr, max_iter, tol);
  const std::size_t J = s.lo.size();
  std::vector<double> cs(J + 1, 0.0);                     // cs[j] = sum_{l >= j} p_l (0-based)
  for (std::size_t j = J; j-- > 0;) cs[j] = cs[j + 1] + (r.p.empty() ? kNaN : r.p[j]);
  auto xlx = [](double v) { return v > 0.0 ? v * std::log(v) : 0.0; };
  std::vector<double> sc(static_cast<std::size_t>(n));
  double mean = 0.0;
  for (int i = 0; i < n; ++i) {
    const double SL = cs[static_cast<std::size_t>(s.a[i])], SR = cs[static_cast<std::size_t>(s.b[i]) + 1];
    sc[i] = (xlx(SL) - xlx(SR)) / (SL - SR);
    mean += sc[i];
  }
  mean /= n;
  out.n.assign(G, 0.0); out.U.assign(G, 0.0); out.sd.assign(G, kNaN);
  std::vector<double> sum(G, 0.0);
  for (int i = 0; i < n; ++i) { out.n[g[i]] += 1.0; sum[g[i]] += sc[i]; }
  for (int k = 0; k < G; ++k) out.U[k] = sum[k] - out.n[k] * mean;
  double vc = 0.0;
  for (int i = 0; i < n; ++i) vc += (sc[i] - mean) * (sc[i] - mean);
  vc /= (n - 1);
  Mat V(static_cast<std::size_t>(G));
  for (int a = 0; a < G; ++a)
    for (int b = 0; b < G; ++b) V(a, b) = vc * ((a == b ? out.n[a] : 0.0) - out.n[a] * out.n[b] / n);
  for (int k = 0; k < G; ++k) out.sd[k] = std::sqrt(V(k, k));
  out.df = G - 1;
  if (G == 2) out.stat = out.U[0] * out.U[0] / V(0, 0);
  else {
    Mat A(static_cast<std::size_t>(G - 1));
    std::vector<double> u(out.U.begin(), out.U.end() - 1);
    for (int a = 0; a < G - 1; ++a) for (int b = 0; b < G - 1; ++b) A(a, b) = V(a, b);
    std::vector<double> x = u;
    if (cholesky(A)) {
      chol_solve(A, x);
      double q = 0.0;
      for (int a = 0; a < G - 1; ++a) q += u[a] * x[a];
      out.stat = q;
    }
  }
  out.p = pchisq_upper(out.stat, out.df);
  return out;
}

// bootstrap of the curve: B resamples of the subjects (seeded, as sc_resample), the NPMLE with the
// resampling counts as case weights on the support of the data, the estimate at the grid times
inline void tb_bootstrap(int n, const double* L, const double* R, std::size_t ng, const double* grid, int B, std::uint64_t seed,
                         int max_iter, double tol, double* out) {
  const TurnbullSupport s = turnbull_support(n, L, R);
  const std::size_t J = s.lo.size();
  std::vector<double> w(static_cast<std::size_t>(n));
  std::uint64_t st = seed;
  for (int b = 0; b < B; ++b) {
    std::fill(w.begin(), w.end(), 0.0);
    for (int i = 0; i < n; ++i) w[draw_index(st, static_cast<std::size_t>(n))] += 1.0;
    const TurnbullResult r = turnbull_fit(s, w.data(), max_iter, tol);
    for (std::size_t k = 0; k < ng; ++k) {
      double e = kNaN, lo, hi;
      if (!r.p.empty()) tb_surv(J, s.lo.data(), s.hi.data(), r.p.data(), grid[k], false, e, lo, hi);
      out[static_cast<std::size_t>(b) * ng + k] = e;
    }
  }
}

}  // namespace survcore
