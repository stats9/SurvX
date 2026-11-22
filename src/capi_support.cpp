// survcore C API: numerical summaries shared by the front ends (survcore/summaries.hpp and the
// modules of the later stages). Same conventions as capi.cpp: 0 on success, -1 bad pointer or
// dimension, -2 non-finite or invalid input, -9 internal error; nothing is thrown across the API.
#include "survcore/survcore.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "survcore/adjust.hpp"
#include "survcore/expand.hpp"
#include "survcore/interval_summary.hpp"
#include "survcore/phtest.hpp"
#include "survcore/timeeffects.hpp"
#include "survcore/toolkit.hpp"
#include "survcore/scaling.hpp"
#include "survcore/summaries.hpp"

using namespace survcore;

namespace {
bool sorted_finite(int m, const double* v) {
  for (int j = 0; j < m; ++j) { if (!std::isfinite(v[j])) return false; if (j && v[j] < v[j - 1]) return false; }
  return true;
}
}  // namespace

extern "C" {

double sc_qnorm(double p) try { return qnorm_p(p); } catch (...) { return kNaN; }
double sc_pnorm(double x) try { return pnorm_std(x); } catch (...) { return kNaN; }
double sc_pchisq_upper(double x, double df) try { return pchisq_upper(x, df); } catch (...) { return kNaN; }
double sc_qgamma(double p, double shape, double rate) try { return qgamma_rate(p, shape, rate); } catch (...) { return kNaN; }

int sc_quantile(int n, const double* x, int np, const double* p, int type, double* out) try {
  if (n < 0 || np < 0 || (n > 0 && !x) || (np > 0 && (!p || !out)) || (type != 1 && type != 7)) return -1;
  std::vector<double> v;
  v.reserve(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) if (!std::isnan(x[i])) v.push_back(x[i]);   // NA dropped, infinite values kept (as R)
  std::sort(v.begin(), v.end());
  for (int j = 0; j < np; ++j) {
    if (!(p[j] >= 0.0 && p[j] <= 1.0)) { out[j] = kNaN; continue; }
    out[j] = type == 1 ? quantile1_sorted(v, p[j]) : quantile7_sorted(v, p[j]);
  }
  return 0;
} catch (...) { return -9; }

int sc_cif_interval(int n, const double* F, const double* se, double zq, double* lower, double* upper) try {
  if (n < 0 || (n > 0 && (!F || !se || !lower || !upper))) return -1;
  for (int i = 0; i < n; ++i) cif_interval(F[i], se[i], zq, lower[i], upper[i]);
  return 0;
} catch (...) { return -9; }

int sc_step_at(int m, const double* t, const double* v, int ns, const double* s, double before, double* out) try {
  if (m < 0 || ns < 0 || (m > 0 && (!t || !v)) || (ns > 0 && (!s || !out))) return -1;
  for (int i = 1; i < m; ++i) if (!(t[i] >= t[i - 1])) return -2;
  for (int j = 0; j < ns; ++j) out[j] = std::isnan(s[j]) ? kNaN : step_at(t, v, static_cast<std::size_t>(m), s[j], before);
  return 0;
} catch (...) { return -9; }

int sc_frailty_lrt(double ll, double ll0, int boundary, double* stat, double* p) try {
  if (!stat || !p) return -1;
  frailty_lrt(ll, ll0, boundary != 0, *stat, *p);
  return 0;
} catch (...) { return -9; }

int sc_frailty_clusters(int n, const int* status, const int* cluster, int C, const double* shape, const double* rate,
                        int boundary, double level, double* size, double* events, double* mean, double* lower,
                        double* upper) try {
  if (n < 0 || C < 1 || (n > 0 && (!status || !cluster)) || !size || !events || !mean || !lower || !upper) return -1;
  if (!boundary && (!shape || !rate)) return -1;
  if (!(level > 0.0 && level < 1.0)) return -2;
  for (int c = 0; c < C; ++c) { size[c] = 0.0; events[c] = 0.0; }
  for (int i = 0; i < n; ++i) {
    if (cluster[i] < 0 || cluster[i] >= C) return -2;
    size[cluster[i]] += 1.0;
    events[cluster[i]] += status[i];
  }
  const double a = (1.0 - level) / 2.0;
  for (int c = 0; c < C; ++c) {
    if (boundary) { mean[c] = lower[c] = upper[c] = 1.0; continue; }
    mean[c] = shape[c] / rate[c];
    lower[c] = qgamma_rate(a, shape[c], rate[c]);
    upper[c] = qgamma_rate(1.0 - a, shape[c], rate[c]);
  }
  return 0;
} catch (...) { return -9; }

// Cox frailty (Breslow jumps hazard at times bh_time, any order): H(t | x) = exp(x'beta) sum_{s <= t} dH(s)
int sc_breslow_cumhaz(int nb, const double* bh_time, const double* bh_hazard, int n, int p, const double* X,
                      const double* beta, int nt, const double* times, double* out) try {
  if (nb < 0 || n < 0 || p < 0 || nt < 0 || (nb > 0 && (!bh_time || !bh_hazard)) || (n > 0 && p > 0 && (!X || !beta)) ||
      (nt > 0 && (!times || (n > 0 && !out))))
    return -1;
  std::vector<double> H0(static_cast<std::size_t>(nt), 0.0);
  for (int j = 0; j < nt; ++j)
    for (int k = 0; k < nb; ++k) if (bh_time[k] <= times[j]) H0[j] += bh_hazard[k];
  for (int i = 0; i < n; ++i) {
    double lp = 0.0;
    for (int r = 0; r < p; ++r) lp += X[static_cast<std::size_t>(i) * p + r] * beta[r];
    const double e = std::exp(lp);
    for (int j = 0; j < nt; ++j) out[static_cast<std::size_t>(i) * nt + j] = e * H0[j];
  }
  return 0;
} catch (...) { return -9; }

// survival of a shared gamma frailty model from the cumulative hazards H (n x nt, row-major):
// mode 0 marginal (1 + theta H)^(-1/theta), 1 conditional exp(-H), 2 cluster (rate_i / (rate_i + H))^shape_i
// (shape, rate: n values, the posterior of the cluster of every row)
int sc_frailty_survival(int n, int nt, const double* H, double theta, int mode, const double* shape, const double* rate,
                        double* out) try {
  if (n < 0 || nt < 0 || (n > 0 && nt > 0 && (!H || !out)) || mode < 0 || mode > 2 || (mode == 2 && (!shape || !rate)))
    return -1;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < nt; ++j) {
      const std::size_t k = static_cast<std::size_t>(i) * nt + j;
      out[k] = mode == 1 ? std::exp(-H[k]) : mode == 0 ? std::pow(1.0 + theta * H[k], -1.0 / theta)
                                                       : std::pow(rate[i] / (rate[i] + H[k]), shape[i]);
    }
  return 0;
} catch (...) { return -9; }

// Nelson-Aalen table from the numbers of events d and at risk Y at the event times: cumulative
// hazard, its standard error sqrt(sum d / Y^2) and the limits H exp(-/+ z se / H)
int sc_na_table(int m, const double* d, const double* Y, double zq, double* H, double* se, double* lower,
                double* upper) try {
  if (m < 0 || (m > 0 && (!d || !Y || !H || !se || !lower || !upper))) return -1;
  double h = 0.0, v = 0.0;
  for (int j = 0; j < m; ++j) {
    h += d[j] / Y[j]; v += d[j] / (Y[j] * Y[j]);
    H[j] = h; se[j] = std::sqrt(v);
    lower[j] = h * std::exp(-zq * se[j] / h); upper[j] = h * std::exp(zq * se[j] / h);
  }
  return 0;
} catch (...) { return -9; }


// Kaplan-Meier limits on the log scale: S / w and min(1, S w), w = exp(z se / S); NaN where S = 0
int sc_km_limits(int m, const double* S, const double* se, double zq, double* lower, double* upper) try {
  if (m < 0 || (m > 0 && (!S || !se || !lower || !upper))) return -1;
  for (int j = 0; j < m; ++j) {
    if (!(S[j] > 0.0)) { lower[j] = upper[j] = kNaN; continue; }
    const double w = std::exp(zq * se[j] / S[j]);
    lower[j] = S[j] / w; upper[j] = std::fmin(1.0, S[j] * w);
  }
  return 0;
} catch (...) { return -9; }

// first time at which v (on the times t) is at or below p: the median of a survival curve and of
// its limits; NaN if never
int sc_first_below(int m, const double* t, const double* v, double p, double* out) try {
  if (m < 0 || !out || (m > 0 && (!t || !v))) return -1;
  *out = kNaN;
  for (int j = 0; j < m; ++j) if (v[j] <= p) { *out = t[j]; break; }
  return 0;
} catch (...) { return -9; }

// numbers of events of cause k and at risk at the distinct event times of cause k (risk set
// {entry < s <= time} with delayed entry; entry may be null); out arrays need room for n; *m the count
int sc_event_counts(int n, const double* time, const int* status, int k, const double* entry, double* et, double* d,
                    double* Y, int* m) try {
  if (n < 0 || !m || (n > 0 && (!time || !status || !et || !d || !Y))) return -1;
  std::vector<double> ev, ts(time, time + n), es;
  for (int i = 0; i < n; ++i) if (status[i] == k) ev.push_back(time[i]);
  std::sort(ev.begin(), ev.end());
  std::sort(ts.begin(), ts.end());
  if (entry) { es.assign(entry, entry + n); std::sort(es.begin(), es.end()); }
  // at risk at s: #{time >= s} - #{entry >= s} (entry < time, so entry >= s implies time >= s)
  auto count_ge = [](const std::vector<double>& v, double s) {
    return static_cast<double>(v.end() - std::lower_bound(v.begin(), v.end(), s));
  };
  int J = 0;
  for (std::size_t a = 0; a < ev.size(); ) {
    std::size_t b = a;
    while (b < ev.size() && ev[b] == ev[a]) ++b;
    et[J] = ev[a]; d[J] = static_cast<double>(b - a);
    Y[J] = count_ge(ts, ev[a]) - (entry ? count_ge(es, ev[a]) : 0.0);
    ++J; a = b;
  }
  *m = J;
  return 0;
} catch (...) { return -9; }


// ---------------- adjusted risks (survcore/adjust.hpp) ----------------
int sc_mlogit(int n, int p, const double* X, const int* a, int K, double* P, int* iterations) try {
  if (n < 1 || p < 1 || K < 2 || !X || !a || !P) return -1;
  for (int i = 0; i < n; ++i) if (a[i] < 0 || a[i] >= K) return -2;
  if (!std::all_of(X, X + static_cast<std::size_t>(n) * p, [](double v) { return std::isfinite(v); })) return -2;
  std::vector<double> Pv;
  int it = 0;
  mlogit_fit(static_cast<std::size_t>(n), static_cast<std::size_t>(p), static_cast<std::size_t>(K), X, a, Pv, it);
  std::copy(Pv.begin(), Pv.end(), P);
  if (iterations) *iterations = it;
  return 0;
} catch (...) { return -9; }

int sc_ipw_weights(int n, int K, const double* P, const int* a, double trim, double* w) try {
  if (n < 1 || K < 2 || !P || !a || !w) return -1;
  if (!(trim >= 0.0 && trim < 0.5)) return -2;
  for (int i = 0; i < n; ++i) if (a[i] < 0 || a[i] >= K) return -2;
  std::vector<double> wv;
  ipw_weights(static_cast<std::size_t>(n), static_cast<std::size_t>(K), P, a, trim, wv);
  std::copy(wv.begin(), wv.end(), w);
  return 0;
} catch (...) { return -9; }

// standardised mean differences of the q columns of Xc between every level l = 1..K-1 and level 0,
// before and after weighting ((K - 1) x q, row-major); NaN where the pooled standard deviation is
// not finite or zero. Weight diagnostics per level: size, effective sample size sum(w)^2 / sum(w^2),
// largest weight, smallest and largest probability P[, l] over all subjects.
int sc_balance(int n, int q, const double* Xc, const int* a, int K, const double* w, const double* P, double* before,
               double* after, double* size, double* ess, double* max_weight, double* ps_min, double* ps_max) try {
  if (n < 1 || q < 0 || K < 2 || !a || !w || !P || !size || !ess || !max_weight || !ps_min || !ps_max) return -1;
  if (q > 0 && (!Xc || !before || !after)) return -1;
  for (int i = 0; i < n; ++i) if (a[i] < 0 || a[i] >= K) return -2;
  auto stats = [&](int lev, int c, double& mean, double& var, double& wmean) {
    double s = 0.0, sw = 0.0, swx = 0.0; int m = 0;
    for (int i = 0; i < n; ++i) if (a[i] == lev) { const double x = Xc[static_cast<std::size_t>(i) * q + c]; s += x; sw += w[i]; swx += x * w[i]; ++m; }
    mean = s / m; wmean = swx / sw;
    double ss = 0.0;
    for (int i = 0; i < n; ++i) if (a[i] == lev) { const double d = Xc[static_cast<std::size_t>(i) * q + c] - mean; ss += d * d; }
    var = m > 1 ? ss / (m - 1) : kNaN;
  };
  for (int l = 1; l < K; ++l)
    for (int c = 0; c < q; ++c) {
      double m0, v0, w0, m1, v1, w1;
      stats(0, c, m0, v0, w0); stats(l, c, m1, v1, w1);
      double sd = std::sqrt((v0 + v1) / 2.0);
      if (!std::isfinite(sd) || sd == 0.0) sd = kNaN;
      before[static_cast<std::size_t>(l - 1) * q + c] = (m1 - m0) / sd;
      after[static_cast<std::size_t>(l - 1) * q + c] = (w1 - w0) / sd;
    }
  for (int l = 0; l < K; ++l) {
    double sw = 0.0, sw2 = 0.0, mx = -INFINITY, cnt = 0.0, pmn = INFINITY, pmx = -INFINITY;
    for (int i = 0; i < n; ++i) {
      const double pl = P[static_cast<std::size_t>(i) * K + l];
      pmn = std::fmin(pmn, pl); pmx = std::fmax(pmx, pl);
      if (a[i] != l) continue;
      sw += w[i]; sw2 += w[i] * w[i]; mx = std::fmax(mx, w[i]); cnt += 1.0;
    }
    size[l] = cnt; ess[l] = sw * sw / sw2; max_weight[l] = mx; ps_min[l] = pmn; ps_max[l] = pmx;
  }
  return 0;
} catch (...) { return -9; }

int sc_weighted_cif(int n, const double* time, const int* status, const double* entry, const double* w, int k, int nt,
                    const double* tt, double* out) try {
  if (n < 0 || nt < 0 || (n > 0 && (!time || !status || !w)) || (nt > 0 && (!tt || !out))) return -1;
  weighted_cif(static_cast<std::size_t>(n), time, status, entry, w, k, static_cast<std::size_t>(nt), tt, out);
  return 0;
} catch (...) { return -9; }

// summaries of nb bootstrap draws of m quantities (draws: nb x m, row-major); see boot_column
int sc_boot_summary(int nb, int m, const double* draws, int finite_only, double level, double* sd, double* lower,
                    double* upper) try {
  if (nb < 0 || m < 0 || (nb > 0 && m > 0 && !draws) || (m > 0 && (!sd || !lower || !upper))) return -1;
  if (!(level > 0.0 && level < 1.0)) return -2;
  const double a = (1.0 - level) / 2.0;
  for (int j = 0; j < m; ++j)
    boot_column(static_cast<std::size_t>(nb), draws + j, static_cast<std::size_t>(m), finite_only != 0, a, sd[j], lower[j], upper[j]);
  return 0;
} catch (...) { return -9; }

int sc_resample(unsigned long long seed, int n, int B, int* out) try {
  if (n < 1 || B < 0 || (B > 0 && !out)) return -1;
  resample_indices(seed, static_cast<std::size_t>(n), static_cast<std::size_t>(B), out);
  return 0;
} catch (...) { return -9; }

int sc_stratified_folds(int n, const int* strata, int folds, unsigned long long seed, int* out) try {
  if (n < 1 || folds < 1 || !strata || !out) return -1;
  stratified_folds(static_cast<std::size_t>(n), strata, folds, seed, out);
  return 0;
} catch (...) { return -9; }


// ---------------- Turnbull NPMLE summaries (survcore/interval_summary.hpp) ----------------
int sc_tb_surv(int J, const double* lo, const double* hi, const double* p, int nt, const double* t, int left, double* est,
               double* lower, double* upper) try {
  if (J < 0 || nt < 0 || (J > 0 && (!lo || !hi || !p)) || (nt > 0 && (!t || !est || !lower || !upper))) return -1;
  for (int k = 0; k < nt; ++k) tb_surv(static_cast<std::size_t>(J), lo, hi, p, t[k], left != 0, est[k], lower[k], upper[k]);
  return 0;
} catch (...) { return -9; }

// time and surv need room for 2 (2 J + 1) values; *m receives the number used
int sc_tb_steps(int J, const double* lo, const double* hi, const double* p, double* time, double* surv, int* m) try {
  if (J < 0 || !time || !surv || !m || (J > 0 && (!lo || !hi || !p))) return -1;
  std::vector<double> t, s;
  tb_steps(static_cast<std::size_t>(J), lo, hi, p, t, s);
  std::copy(t.begin(), t.end(), time); std::copy(s.begin(), s.end(), surv);
  *m = static_cast<int>(t.size());
  return 0;
} catch (...) { return -9; }

int sc_tb_quantile(int m, const double* time, const double* surv, double p, double* out) try {
  if (m < 0 || !out || (m > 0 && (!time || !surv))) return -1;
  *out = tb_quantile(static_cast<std::size_t>(m), time, surv, p);
  return 0;
} catch (...) { return -9; }

int sc_tb_test(int n, const double* L, const double* R, const int* g, int G, int max_iter, double tol, double* stat, int* df,
               double* p, double* size, double* U, double* sd) try {
  if (n < 2 || G < 2 || !L || !R || !g || !stat || !df || !p || !size || !U || !sd) return -1;
  for (int i = 0; i < n; ++i) {
    if (g[i] < 0 || g[i] >= G) return -2;
    if (!std::isfinite(L[i]) || L[i] < 0.0 || std::isnan(R[i]) || R[i] < L[i]) return -2;
  }
  const TurnbullTest r = tb_test(n, L, R, g, G, max_iter, tol);
  *stat = r.stat; *df = r.df; *p = r.p;
  for (int k = 0; k < G; ++k) { size[k] = r.n[k]; U[k] = r.U[k]; sd[k] = r.sd[k]; }
  return 0;
} catch (...) { return -9; }

int sc_tb_bootstrap(int n, const double* L, const double* R, int ng, const double* grid, int B, unsigned long long seed,
                    int max_iter, double tol, double* out) try {
  if (n < 1 || ng < 0 || B < 0 || !L || !R || (ng > 0 && !grid) || (B > 0 && ng > 0 && !out)) return -1;
  for (int i = 0; i < n; ++i) if (!std::isfinite(L[i]) || L[i] < 0.0 || std::isnan(R[i]) || R[i] < L[i]) return -2;
  tb_bootstrap(n, L, R, static_cast<std::size_t>(ng), grid, B, seed, max_iter, tol, out);
  return 0;
} catch (...) { return -9; }

// percentile limits at `times` of bootstrap curves stored on `grid` (draws B x ng): the column of
// a grid time within 1e-8 max(1, t), otherwise the curves interpolated linearly between the
// neighbouring grid times (NaN beyond the grid); limits from at least min_n finite values
int sc_grid_boot_limits(int ng, const double* grid, int B, const double* draws, int nt, const double* times, double level,
                        int min_n, double* lower, double* upper) try {
  if (ng < 1 || B < 0 || nt < 0 || !grid || (B > 0 && !draws) || (nt > 0 && (!times || !lower || !upper))) return -1;
  const double a = (1.0 - level) / 2.0;
  const double gmax = *std::max_element(grid, grid + ng);
  for (int m = 0; m < nt; ++m) {
    lower[m] = upper[m] = kNaN;
    const double t = times[m];
    int j = 0;
    for (int k = 1; k < ng; ++k) if (std::fabs(grid[k] - t) < std::fabs(grid[j] - t)) j = k;
    std::vector<double> v;
    if (std::fabs(grid[j] - t) < 1e-8 * std::fmax(1.0, t)) {
      for (int b = 0; b < B; ++b) v.push_back(draws[static_cast<std::size_t>(b) * ng + j]);
    } else {
      int j1 = -1, j2 = -1;
      for (int k = 0; k < ng; ++k) { if (grid[k] < t) j1 = k; if (grid[k] > t && j2 < 0) j2 = k; }
      if (j1 < 0) j1 = 0;
      if (j2 < 0) j2 = ng - 1;
      if (j1 == j2 || !std::isfinite(t) || t > gmax) continue;
      const double f = (t - grid[j1]) / (grid[j2] - grid[j1]);
      for (int b = 0; b < B; ++b)
        v.push_back((1.0 - f) * draws[static_cast<std::size_t>(b) * ng + j1] + f * draws[static_cast<std::size_t>(b) * ng + j2]);
    }
    std::vector<double> q;
    for (double y : v) if (std::isfinite(y)) q.push_back(y);
    if (static_cast<int>(q.size()) < min_n) continue;
    std::sort(q.begin(), q.end());
    lower[m] = quantile7_sorted(q, a); upper[m] = quantile7_sorted(q, 1.0 - a);
  }
  return 0;
} catch (...) { return -9; }

// percentile limits of every column of B x m draws from its finite values (NaN with fewer than min_n)
int sc_column_limits(int B, int m, const double* draws, double level, int min_n, double* lower, double* upper) try {
  if (B < 0 || m < 0 || (B > 0 && m > 0 && !draws) || (m > 0 && (!lower || !upper))) return -1;
  const double a = (1.0 - level) / 2.0;
  for (int j = 0; j < m; ++j) {
    std::vector<double> q;
    for (int b = 0; b < B; ++b) { const double y = draws[static_cast<std::size_t>(b) * m + j]; if (std::isfinite(y)) q.push_back(y); }
    lower[j] = upper[j] = kNaN;
    if (static_cast<int>(q.size()) < min_n) continue;
    std::sort(q.begin(), q.end());
    lower[j] = quantile7_sorted(q, a); upper[j] = quantile7_sorted(q, 1.0 - a);
  }
  return 0;
} catch (...) { return -9; }

// interval for the median of bootstrap curves on the grid (draws B x ng): the median of every curve
// by tb_quantile, limits when at least half of the B curves and at least 10 reach 0.5
int sc_tb_median_limits(int ng, const double* grid, int B, const double* draws, double level, double* lower, double* upper) try {
  if (ng < 1 || B < 0 || !grid || (B > 0 && !draws) || !lower || !upper) return -1;
  std::vector<double> q;
  for (int b = 0; b < B; ++b) {
    const double md = tb_quantile(static_cast<std::size_t>(ng), grid, draws + static_cast<std::size_t>(b) * ng, 0.5);
    if (std::isfinite(md)) q.push_back(md);
  }
  *lower = *upper = kNaN;
  if (static_cast<double>(q.size()) >= 0.5 * B && q.size() >= 10) {
    std::sort(q.begin(), q.end());
    const double a = (1.0 - level) / 2.0;
    *lower = quantile7_sorted(q, a); *upper = quantile7_sorted(q, 1.0 - a);
  }
  return 0;
} catch (...) { return -9; }

// every row of Y (nr x nx) interpolated linearly at xout from its finite values, as R's
// approx(rule = 1); a row with fewer than two finite values gives NaN
int sc_approx_rows(int nr, int nx, const double* x, const double* Y, int nout, const double* xout, double* out) try {
  if (nr < 0 || nx < 0 || nout < 0 || (nx > 0 && !x) || (nr > 0 && nx > 0 && !Y) || (nout > 0 && (!xout || (nr > 0 && !out)))) return -1;
  for (int r = 0; r < nr; ++r) {
    std::vector<double> xs, ys;
    for (int k = 0; k < nx; ++k) { const double y = Y[static_cast<std::size_t>(r) * nx + k]; if (std::isfinite(y)) { xs.push_back(x[k]); ys.push_back(y); } }
    for (int k = 0; k < nout; ++k)
      out[static_cast<std::size_t>(r) * nout + k] = xs.size() < 2 ? kNaN : approx_linear(xs, ys, xout[k]);
  }
  return 0;
} catch (...) { return -9; }


// ---------------- prediction measures and model comparison helpers ----------------
// Brier score (and that of the model without covariates), IPA, time-dependent AUC and calibration
// error of the predicted risks (n x nt, row-major) of cause `interest` at every time; NaN where a
// measure cannot be computed at a time
int sc_measures(int n, const double* time, const int* cause, int interest, int nt, const double* t, const double* risk,
                int groups, double* brier, double* brier_null, double* ipa, double* auc, double* calib) try {
  if (n < 1 || nt < 0 || !time || !cause || (nt > 0 && (!t || !risk || !brier || !brier_null || !ipa || !auc || !calib))) return -1;
  std::vector<double> col(static_cast<std::size_t>(n)), mp(static_cast<std::size_t>(std::max(groups, 1))), ob(mp.size()), sz(mp.size());
  for (int j = 0; j < nt; ++j) {
    for (int i = 0; i < n; ++i) col[i] = risk[static_cast<std::size_t>(i) * nt + j];
    double b, b0, ip, a, e; int m = 0;
    if (sc_brier(n, time, cause, interest, col.data(), t[j], &b, &b0, &ip) == 0) { brier[j] = b; brier_null[j] = b0; ipa[j] = ip; }
    else brier[j] = brier_null[j] = ipa[j] = kNaN;
    auc[j] = sc_auc(n, time, cause, interest, col.data(), t[j], &a) == 0 ? a : kNaN;
    calib[j] = sc_calibration(n, time, cause, interest, col.data(), t[j], groups, mp.data(), ob.data(), sz.data(), &e, &m) == 0 ? e : kNaN;
  }
  return 0;
} catch (...) { return -9; }

// trapezoidal integral of y over (0, grid] with y(0) = 0, divided by the last grid point; NaN if a
// value is not finite
int sc_integrate_grid(int m, const double* grid, const double* y, double* out) try {
  if (m < 1 || !grid || !y || !out) return -1;
  double s = 0.0, x0 = 0.0, v0 = 0.0, mx = grid[0];
  for (int j = 0; j < m; ++j) if (!std::isfinite(y[j])) { *out = kNaN; return 0; }
  for (int j = 0; j < m; ++j) { s += (grid[j] - x0) * (y[j] + v0) / 2.0; x0 = grid[j]; v0 = y[j]; mx = std::fmax(mx, grid[j]); }
  *out = s / mx;
  return 0;
} catch (...) { return -9; }

// non-parametric probability of the event (one event type: 1 - Kaplan-Meier) or cumulative
// incidence of `cause` (Aalen-Johansen) at the times `at`
int sc_np_curve(int n, const double* time, const int* status, int cause, int nt, const double* at, double* out) try {
  if (n < 1 || nt < 0 || !time || !status || (nt > 0 && (!at || !out))) return -1;
  int K = 0;
  for (int i = 0; i < n; ++i) K = std::max(K, status[i]);
  std::vector<double> ot(static_cast<std::size_t>(n)), sv(static_cast<std::size_t>(n));
  int m = 0;
  if (K <= 1) {
    std::vector<double> nr(static_cast<std::size_t>(n)), ne(static_cast<std::size_t>(n)), se(static_cast<std::size_t>(n));
    if (int rc = sc_kaplan_meier(n, time, status, nullptr, ot.data(), nr.data(), ne.data(), sv.data(), se.data(), &m)) return rc;
    for (int j = 0; j < nt; ++j) out[j] = 1.0 - step_at(ot.data(), sv.data(), static_cast<std::size_t>(m), at[j], 1.0);
    return 0;
  }
  if (cause < 1 || cause > K) return -1;
  std::vector<double> cif(static_cast<std::size_t>(n) * K), cse(cif.size());
  if (int rc = sc_aalen_johansen(n, time, status, K, ot.data(), sv.data(), cif.data(), cse.data(), &m)) return rc;
  std::vector<double> ck(static_cast<std::size_t>(m));
  for (int j = 0; j < m; ++j) ck[j] = cif[static_cast<std::size_t>(j) * K + (cause - 1)];
  for (int j = 0; j < nt; ++j) out[j] = step_at(ot.data(), ck.data(), static_cast<std::size_t>(m), at[j], 0.0);
  return 0;
} catch (...) { return -9; }


// ranking of M models on a cross-validated loss (lower is better) with the one-standard-error rule
// (Hastie et al. 2009): fold_loss M x folds (row-major). For the models with a finite loss: the
// difference from the best, the paired standard error of that difference over the folds with a
// finite difference (NaN with fewer than two), equivalent = within one standard error, rank 1 for
// the equivalent models and 1 + (number equivalent), ... for the others in order of the loss;
// rank 0 (and NaN) for models without a finite loss.
int sc_rank_one_se(int M, int folds, const double* loss, const double* fold_loss, int* rank, double* diff, double* se,
                   int* equivalent) try {
  if (M < 1 || folds < 1 || !loss || !fold_loss || !rank || !diff || !se || !equivalent) return -1;
  int best = -1;
  for (int m = 0; m < M; ++m) {
    rank[m] = 0; diff[m] = se[m] = kNaN; equivalent[m] = 0;
    if (std::isfinite(loss[m]) && (best < 0 || loss[m] < loss[best])) best = m;
  }
  if (best < 0) return 0;
  std::vector<int> rest;
  int neq = 0;
  for (int m = 0; m < M; ++m) {
    if (!std::isfinite(loss[m])) continue;
    std::vector<double> dk;
    for (int k = 0; k < folds; ++k) {
      const double d = fold_loss[static_cast<std::size_t>(m) * folds + k] - fold_loss[static_cast<std::size_t>(best) * folds + k];
      if (std::isfinite(d)) dk.push_back(d);
    }
    diff[m] = loss[m] - loss[best];
    if (m == best) se[m] = 0.0;
    else if (dk.size() >= 2) {
      double mu = 0.0;
      for (double d : dk) mu += d;
      mu /= static_cast<double>(dk.size());
      double ss = 0.0;
      for (double d : dk) ss += (d - mu) * (d - mu);
      se[m] = std::sqrt(ss / static_cast<double>(dk.size() - 1)) / std::sqrt(static_cast<double>(dk.size()));
    }
    equivalent[m] = m == best || (std::isfinite(se[m]) && diff[m] <= se[m]);
    if (equivalent[m]) { rank[m] = 1; ++neq; } else rest.push_back(m);
  }
  std::stable_sort(rest.begin(), rest.end(), [&](int a, int b) { return loss[a] < loss[b]; });
  for (std::size_t q = 0; q < rest.size(); ++q) rank[rest[q]] = neq + 1 + static_cast<int>(q);
  return 0;
} catch (...) { return -9; }


// ---------------- test of proportional hazards (survcore/phtest.hpp) ----------------
// status 0/1; w, entry, strata may be null; transform: 0 identity, 1 log, 2 rank, 3 Kaplan-Meier,
// 4 weighted Kaplan-Meier; term_of: term (0..nterm-1) of every column. Outputs: chi-square and p
// per column (p) and per term (nterm, with df), global; for the nd events (*nd, by stratum and time):
// row (0-based), time, transform g and scaled Schoenfeld residuals (nd x p). Row outputs need room
// for n events. Returns -4 when the risk-set information is singular or there is no event.
int sc_ph_test(int n, int p, const double* time, const int* status, const double* X, const double* beta, const double* vcov,
               const double* w, const double* entry, const int* strata, int transform, const int* term_of, int nterm,
               double* chi_col, double* p_col, double* chi_term, double* df_term, double* p_term, double* chi_glob,
               double* p_glob, int* nd, int* rows, double* ev_time, double* g, double* scaled) try {
  if (n < 1 || p < 1 || nterm < 1 || !time || !status || !X || !beta || !vcov || !term_of || !chi_col || !p_col || !chi_term ||
      !df_term || !p_term || !chi_glob || !p_glob || !nd || !rows || !ev_time || !g || !scaled || transform < 0 || transform > 4)
    return -1;
  for (int a = 0; a < p; ++a) if (term_of[a] < 0 || term_of[a] >= nterm) return -1;
  const PhTest r = ph_test(static_cast<std::size_t>(n), static_cast<std::size_t>(p), time, status, X, beta, vcov, w, entry, strata,
                           transform, term_of, static_cast<std::size_t>(nterm));
  if (!r.ok) return -4;
  for (int a = 0; a < p; ++a) { chi_col[a] = r.chi_col[a]; p_col[a] = pchisq_upper(r.chi_col[a], 1.0); }
  for (int t = 0; t < nterm; ++t) { chi_term[t] = r.chi_term[t]; df_term[t] = r.df_term[t]; p_term[t] = pchisq_upper(r.chi_term[t], r.df_term[t]); }
  *chi_glob = r.chi_glob; *p_glob = pchisq_upper(r.chi_glob, static_cast<double>(p));
  *nd = static_cast<int>(r.rows.size());
  for (std::size_t q = 0; q < r.rows.size(); ++q) {
    rows[q] = static_cast<int>(r.rows[q]); ev_time[q] = r.time[q]; g[q] = r.g0[q];
    for (int a = 0; a < p; ++a) scaled[q * static_cast<std::size_t>(p) + a] = r.scaled[q * static_cast<std::size_t>(p) + a];
  }
  return 0;
} catch (...) { return -9; }


// Holm's step-down adjustment of p-values (NaN entries are left out and stay NaN), as R's p.adjust
int sc_p_adjust_holm(int n, const double* p, double* out) try {
  if (n < 0 || (n > 0 && (!p || !out))) return -1;
  std::vector<int> idx;
  for (int i = 0; i < n; ++i) { out[i] = kNaN; if (!std::isnan(p[i])) idx.push_back(i); }
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return p[a] < p[b]; });
  const double m = static_cast<double>(idx.size());
  double run = -INFINITY;
  for (std::size_t q = 0; q < idx.size(); ++q) {
    run = std::fmax(run, (m - static_cast<double>(q)) * p[idx[q]]);
    out[idx[q]] = std::fmin(1.0, run);
  }
  return 0;
} catch (...) { return -9; }


// ---------------- time-varying effects (survcore/timeeffects.hpp) ----------------

int sc_te_split_size(int n, const double* start, const double* stop, int m, const double* et, double* out) try {
  if (n < 0 || m < 0 || !out || (n > 0 && (!start || !stop)) || (m > 0 && !et) || !sorted_finite(m, et)) return -1;
  *out = split_size(static_cast<std::size_t>(n), start, stop, std::vector<double>(et, et + m));
  return 0;
} catch (...) { return -9; }

// cap: room in the outputs (at least sc_te_split_size); *nrow receives the number of rows; idx 0-based
int sc_te_split(int n, const double* start, const double* stop, const int* status, int m, const double* et, int cap, int* idx,
                double* s0, double* s1, int* st, int* nrow) try {
  if (n < 0 || m < 0 || cap < 0 || !nrow || (n > 0 && (!start || !stop || !status)) || (m > 0 && !et) || !sorted_finite(m, et)) return -1;
  std::vector<int> I, S; std::vector<double> A, Bv;
  split_rows(static_cast<std::size_t>(n), start, stop, status, std::vector<double>(et, et + m), I, A, Bv, S);
  if (static_cast<int>(I.size()) > cap || (!I.empty() && (!idx || !s0 || !s1 || !st))) return -1;
  for (std::size_t r = 0; r < I.size(); ++r) { idx[r] = I[r]; s0[r] = A[r]; s1[r] = Bv[r]; st[r] = S[r]; }
  *nrow = static_cast<int>(I.size());
  return 0;
} catch (...) { return -9; }

// out: nt x q (row-major), q = 1 (log, linear), nbreak (steps) or nknot + 1 (spline)
int sc_te_basis(int form, int nt, const double* t, double tmin, double tm, int nbreak, const double* breaks, int nknot,
                const double* knots, double bl, double br, double* out) try {
  if (form < 1 || form > 4 || nt < 0 || nbreak < 0 || nknot < 0 || (nt > 0 && (!t || !out)) || (nbreak > 0 && !breaks) ||
      (nknot > 0 && !knots))
    return -1;
  if (form == TEF_STEP && nbreak < 1) return -1;
  if (form == TEF_SPLINE && !(bl < br)) return -2;
  te_basis(form, static_cast<std::size_t>(nt), t, tmin, tm, breaks, static_cast<std::size_t>(nbreak), knots,
           static_cast<std::size_t>(nknot), bl, br, out);
  return 0;
} catch (...) { return -9; }

int sc_te_effect(int nt, int q, const double* B, int P, const double* b, const double* V, const int* ix, double zq, double* est,
                 double* se, double* hr, double* lower, double* upper, double* p) try {
  if (nt < 0 || q < 0 || P < 1 || !b || !V || !ix || (nt > 0 && (!est || !se || !hr || !lower || !upper || !p)) || (nt > 0 && q > 0 && !B))
    return -1;
  for (int c = 0; c <= q; ++c) if ((c > 0 || ix[c] >= 0) && (ix[c] < 0 || ix[c] >= P)) return -1;
  if (ix[0] < 0 && q < 1) return -1;
  te_effect(static_cast<std::size_t>(nt), static_cast<std::size_t>(q), B, static_cast<std::size_t>(P), b, V, ix, zq, est, se, hr,
            lower, upper, p);
  return 0;
} catch (...) { return -9; }

// dH: K x nr x m (cause, row, event time), row-major; out nr x nt
int sc_cum_predict(int mode, int nr, int m, const double* tj, int K, const double* dH, int cause, int product_limit, int nt,
                   const double* times, double* out) try {
  if (mode < 0 || mode > 2 || nr < 0 || m < 0 || K < 1 || nt < 0 || (m > 0 && (!tj || (nr > 0 && !dH))) || (nt > 0 && (!times || (nr > 0 && !out))))
    return -1;
  if (mode == 2 && (cause < 0 || cause >= K)) return -1;
  if (!sorted_finite(m, tj)) return -2;
  cum_predict(mode, static_cast<std::size_t>(nr), static_cast<std::size_t>(m), tj, static_cast<std::size_t>(K), dH, cause,
              product_limit != 0, static_cast<std::size_t>(nt), times, out);
  return 0;
} catch (...) { return -9; }


// ---- data expansions ---------------------------------------------------------------------------
// Fine-Gray data with delayed entry: first call with cap = 0 to get *nrow, then with room for it
int sc_fg_expand(int n, const double* time, const int* status, const double* entry, int cause, int cap, int* id,
                 double* start, double* stop, int* st, double* w, int* nrow) try {
  if (n < 0 || cap < 0 || !nrow || cause < 1 || (n > 0 && (!time || !status))) return -1;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(time[i]) || status[i] < 0) return -1;
    if (entry && (!std::isfinite(entry[i]) || entry[i] < 0 || entry[i] >= time[i])) return -1;
  }
  const FGRows r = fg_expand(static_cast<std::size_t>(n), time, status, entry, cause);
  *nrow = static_cast<int>(r.id.size());
  if (cap == 0) return 0;
  if (cap < *nrow || !id || !start || !stop || !st || !w) return -1;
  for (std::size_t k = 0; k < r.id.size(); ++k) { id[k] = r.id[k]; start[k] = r.start[k]; stop[k] = r.stop[k]; st[k] = r.status[k]; w[k] = r.w[k]; }
  return 0;
} catch (...) { return -9; }

// merge times closer than tol (in place); *changed = 1 when any were merged
int sc_timefix(int n, double* v, double tol, int* changed) try {
  if (n < 0 || !changed || (n > 0 && !v) || !(tol >= 0)) return -1;
  *changed = time_fix(static_cast<std::size_t>(n), v, tol) ? 1 : 0;
  return 0;
} catch (...) { return -9; }

// robust variance V (sum over clusters of U_c U_c') V; cl: 0-based codes < nc, or null
int sc_sandwich(int n, int p, const double* U, const int* cl, int nc, const double* V, double* out) try {
  if (n < 0 || p < 1 || !V || !out || (n > 0 && !U) || (cl && nc < 1)) return -1;
  if (cl) for (int i = 0; i < n; ++i) if (cl[i] < 0 || cl[i] >= nc) return -1;
  sandwich(static_cast<std::size_t>(n), static_cast<std::size_t>(p), U, cl, static_cast<std::size_t>(cl ? nc : 0), V, out);
  return 0;
} catch (...) { return -9; }


// ---- source of random indices ------------------------------------------------------------------
void sc_set_index_source(sc_index_fn fn) { index_source() = fn; }


// covariance of the complete rows of nb x p draws (row-major); *nok: rows used
int sc_cov_complete(int nb, int p, const double* X, double* out, int* nok) try {
  if (nb < 0 || p < 1 || !out || !nok || (nb > 0 && !X)) return -1;
  std::size_t k = 0;
  cov_complete(static_cast<std::size_t>(nb), static_cast<std::size_t>(p), X, out, k);
  *nok = static_cast<int>(k);
  return 0;
} catch (...) { return -9; }


// smoother of the residual plots (timeeffects.hpp ns_smooth); -4 when the fit is not possible
int sc_ns_smooth(int n, const double* x, const double* y, int df, int ng, const double* grid, double z, double* fit, double* lower,
                 double* upper) try {
  if (n < 0 || ng < 0 || df < 1 || (n > 0 && (!x || !y)) || (ng > 0 && (!grid || !fit || !lower || !upper))) return -1;
  for (int i = 0; i < n; ++i) if (!std::isfinite(x[i]) || !std::isfinite(y[i])) return -1;
  for (int g = 0; g < ng; ++g) if (!std::isfinite(grid[g])) return -1;
  return ns_smooth(static_cast<std::size_t>(n), x, y, df, static_cast<std::size_t>(ng), grid, z, fit, lower, upper) ? 0 : -4;
} catch (...) { return -9; }


// ---- elementwise prediction formulas -------------------------------------------------------------
// Fine-Gray prediction: out (nr x nt, row-major) = 1 - exp(-exp(x_i' beta) H_j)
int sc_subdist_cif(int nr, int p, const double* X, const double* beta, int nt, const double* H, double* out) try {
  if (nr < 0 || p < 1 || nt < 0 || !beta || (nr > 0 && !X) || (nt > 0 && !H) || (nr > 0 && nt > 0 && !out)) return -1;
  for (int i = 0; i < nr; ++i) {
    double lp = 0.0;
    for (int j = 0; j < p; ++j) lp += X[static_cast<std::size_t>(i) * p + j] * beta[j];
    const double r = std::exp(lp);
    for (int j = 0; j < nt; ++j) out[static_cast<std::size_t>(i) * nt + j] = -std::expm1(-r * H[j]);
  }
  return 0;
} catch (...) { return -9; }

// limits on the log scale: est exp(-+ z se / est); NaN where est or se is not finite or est <= 0
int sc_log_limits(int n, const double* est, const double* se, double z, double* lower, double* upper) try {
  if (n < 0 || (n > 0 && (!est || !se || !lower || !upper))) return -1;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(est[i]) || !std::isfinite(se[i]) || !(est[i] > 0)) { lower[i] = upper[i] = std::nan(""); continue; }
    const double e = z * se[i] / est[i];
    lower[i] = est[i] * std::exp(-e); upper[i] = est[i] * std::exp(e);
  }
  return 0;
} catch (...) { return -9; }


// ---- toolkit: the small computations on results (survcore/toolkit.hpp) -------------------------
int sc_elementwise(int code, int n, const double* x, const double* y, double* out) try {
  if (!ew_valid(code) || n < 0 || (n > 0 && (!x || !out)) || (n > 0 && ew_binary(code) && !y)) return -1;
  for (int i = 0; i < n; ++i) out[i] = ew1(code, x[i], ew_binary(code) ? y[i] : 0.0);
  return 0;
} catch (...) { return -9; }

int sc_reduce(int code, int n, const double* x, int na_rm, double* out) try {
  if (code < RD_SUM || code > RD_MIN || n < 0 || !out || (n > 0 && !x)) return -1;
  *out = reduce(code, static_cast<std::size_t>(n), x, 1, na_rm != 0);
  return 0;
} catch (...) { return -9; }

// by_row = 0: means of the nc columns of M (nr x nc, row-major); 1: means of the nr rows
int sc_means(int nr, int nc, const double* M, int by_row, int na_rm, double* out) try {
  if (nr < 0 || nc < 0 || (nr > 0 && nc > 0 && !M) || ((by_row ? nr : nc) > 0 && !out)) return -1;
  if (by_row) for (int i = 0; i < nr; ++i) out[i] = r_colmean(static_cast<std::size_t>(nc), M + static_cast<std::size_t>(i) * nc, 1, na_rm != 0);
  else for (int j = 0; j < nc; ++j) out[j] = r_colmean(static_cast<std::size_t>(nr), M + j, static_cast<std::size_t>(nc), na_rm != 0);
  return 0;
} catch (...) { return -9; }

// means of x within groups g (0-based codes < G), NaN left out (as aggregate(), which drops NA)
int sc_group_means(int n, const int* g, const double* x, int G, double* out) try {
  if (n < 0 || G < 1 || !out || (n > 0 && (!g || !x))) return -1;
  std::vector<std::vector<double>> v(static_cast<std::size_t>(G));
  for (int i = 0; i < n; ++i) { if (g[i] < 0 || g[i] >= G) return -1; if (!std::isnan(x[i])) v[static_cast<std::size_t>(g[i])].push_back(x[i]); }
  for (int k = 0; k < G; ++k) out[k] = r_mean_of(v[static_cast<std::size_t>(k)]);
  return 0;
} catch (...) { return -9; }

int sc_cumsum_rows(int nr, int nc, const double* M, double cap, double* out) try {
  if (nr < 0 || nc < 0 || (nr > 0 && nc > 0 && (!M || !out))) return -1;
  cumsum_rows(static_cast<std::size_t>(nr), static_cast<std::size_t>(nc), M, cap, out);
  return 0;
} catch (...) { return -9; }

// X beta for X nr x p (row-major)
int sc_linpred(int nr, int p, const double* X, const double* beta, double* out) try {
  if (nr < 0 || p < 0 || (nr > 0 && !out) || (nr > 0 && p > 0 && (!X || !beta))) return -1;
  for (int i = 0; i < nr; ++i) {
    double s = 0.0;
    for (int j = 0; j < p; ++j) s += X[static_cast<std::size_t>(i) * p + j] * beta[j];
    out[i] = s;
  }
  return 0;
} catch (...) { return -9; }

// proportional-hazards prediction from a cumulative baseline: mode 0 survival exp(-e^{lp} H),
// 1 cumulative incidence 1 - exp(-e^{lp} H); lp = X beta (X may be null with p = 0); H a vector of
// nt values, or (by_row) a nr x nt matrix with the baseline of each row (strata)
int sc_ph_predict(int mode, int nr, int p, const double* X, const double* beta, int nt, const double* H, int by_row, double* out) try {
  if ((mode != 0 && mode != 1) || nr < 0 || p < 0 || nt < 0 || (nr > 0 && p > 0 && (!X || !beta)) || (nt > 0 && !H) ||
      (nr > 0 && nt > 0 && !out))
    return -1;
  for (int i = 0; i < nr; ++i) {
    double lp = 0.0;
    for (int j = 0; j < p; ++j) lp += X[static_cast<std::size_t>(i) * p + j] * beta[j];
    const double r = std::exp(lp);
    for (int j = 0; j < nt; ++j) {
      const double h = by_row ? H[static_cast<std::size_t>(i) * nt + j] : H[j];
      const double S = std::exp(-r * h);
      out[static_cast<std::size_t>(i) * nt + j] = mode == 0 ? S : 1.0 - S;
    }
  }
  return 0;
} catch (...) { return -9; }

int sc_np_curve_inc(int m, const double* d, const double* Y, int mode, double* out) try {
  if (m < 0 || mode < 0 || mode > 2 || (m > 0 && (!d || !Y || !out))) return -1;
  np_curve_inc(static_cast<std::size_t>(m), d, Y, mode, out);
  return 0;
} catch (...) { return -9; }

// AIC -2 ll + 2 k and BIC -2 ll + log(nobs) k, elementwise
int sc_info_criteria(int n, const double* ll, const double* k, const double* nobs, double* aic, double* bic) try {
  if (n < 0 || (n > 0 && (!ll || !k || !aic || !bic))) return -1;
  for (int i = 0; i < n; ++i) {
    aic[i] = -2.0 * ll[i] + 2.0 * k[i];
    bic[i] = nobs ? -2.0 * ll[i] + std::log(nobs[i]) * k[i] : std::nan("");
  }
  return 0;
} catch (...) { return -9; }

// likelihood ratio tests 2 (ll1 - ll0) on df degrees of freedom, elementwise; p NaN when df < 1
int sc_lr_test(int n, const double* ll1, const double* ll0, const double* df, double* stat, double* p) try {
  if (n < 0 || (n > 0 && (!ll1 || !ll0 || !df || !stat || !p))) return -1;
  for (int i = 0; i < n; ++i) {
    stat[i] = 2.0 * (ll1[i] - ll0[i]);
    p[i] = df[i] >= 1 && std::isfinite(stat[i]) ? pchisq_upper(std::fmax(stat[i], 0.0), df[i]) : std::nan("");
  }
  return 0;
} catch (...) { return -9; }

// Wald table: limits est -+ z se, z statistic, two-sided p, and where expo[i] the exponentiated
// estimate and limits (NaN otherwise)
int sc_wald_table(int n, const double* est, const double* se, double z, const int* expo, double* lo, double* hi, double* zs,
                  double* p, double* eest, double* elo, double* ehi) try {
  if (n < 0 || (n > 0 && (!est || !se || !expo || !lo || !hi || !zs || !p || !eest || !elo || !ehi))) return -1;
  for (int i = 0; i < n; ++i) {
    lo[i] = est[i] - z * se[i]; hi[i] = est[i] + z * se[i];
    zs[i] = est[i] / se[i];
    p[i] = std::isnan(zs[i]) ? zs[i] : 2.0 * pnorm_std(-std::fabs(zs[i]));
    eest[i] = expo[i] ? std::exp(est[i]) : std::nan("");
    elo[i] = expo[i] ? std::exp(lo[i]) : std::nan("");
    ehi[i] = expo[i] ? std::exp(hi[i]) : std::nan("");
  }
  return 0;
} catch (...) { return -9; }

// standard errors sqrt(diag(V)), V p x p; NaN for a negative or missing variance
int sc_se_diag(int p, const double* V, double* out) try {
  if (p < 0 || (p > 0 && (!V || !out))) return -1;
  for (int j = 0; j < p; ++j) { const double v = V[static_cast<std::size_t>(j) * p + j]; out[j] = v >= 0 ? std::sqrt(v) : std::nan(""); }
  return 0;
} catch (...) { return -9; }

// scaled Schoenfeld residuals nd (S V) + beta, S nd x p row-major, V p x p
int sc_scaled_schoenfeld(int nd, int p, const double* S, const double* V, const double* beta, double* out) try {
  if (nd < 0 || p < 1 || !V || !beta || (nd > 0 && (!S || !out))) return -1;
  for (int i = 0; i < nd; ++i)
    for (int k = 0; k < p; ++k) {
      double s = 0.0;
      for (int j = 0; j < p; ++j) s += S[static_cast<std::size_t>(i) * p + j] * V[static_cast<std::size_t>(j) * p + k];
      out[static_cast<std::size_t>(i) * p + k] = s * nd + beta[k];
    }
  return 0;
} catch (...) { return -9; }


// grids: kind 0 seq(from, to, length.out = n); 1 seq(from, to, by); 2 pretty(c(from, to), n).
// First call with cap = 0 for *m, then with room for it.
int sc_grid(int kind, double from, double to, double by, int n, int cap, double* out, int* m) try {
  if (kind < 0 || kind > 2 || cap < 0 || !m || !std::isfinite(from) || !std::isfinite(to)) return -1;
  if ((kind == 0 || kind == 2) && n < 0) return -1;
  if (kind == 1 && !(std::isfinite(by) && by != 0)) return -1;
  const std::vector<double> g = kind == 0 ? seq_length(from, to, static_cast<std::size_t>(n)) : kind == 1 ? seq_by(from, to, by)
                                                                                                      : pretty_r(std::fmin(from, to), std::fmax(from, to), n);
  *m = static_cast<int>(g.size());
  if (cap == 0) return 0;
  if (cap < *m || !out) return -1;
  for (std::size_t k = 0; k < g.size(); ++k) out[k] = g[k];
  return 0;
} catch (...) { return -9; }


// hazard increments of rows on a grid: the baseline (time, cumulative hazard, stratum; sorted by
// time within each stratum) gives increments dH0 at its times; row i of stratum s gets
// dH0 * exp(LP) at the matching grid time (LP: nr x m row-major when lp_by_time, else one value
// per row). Grid times without a jump of the row's stratum stay 0. -1 when a baseline time is not
// on the grid.
int sc_hazard_increments(int nb, const double* btime, const double* bcum, const int* bstr, int nr, const int* rstr, int m,
                         const double* tj, const double* LP, int lp_by_time, double* out) try {
  if (nb < 0 || nr < 0 || m < 0 || (nb > 0 && (!btime || !bcum || !bstr)) || (nr > 0 && (!rstr || !LP)) || (m > 0 && !tj) ||
      (nr > 0 && m > 0 && !out))
    return -1;
  for (std::size_t k = 0; k < static_cast<std::size_t>(nr) * m; ++k) out[k] = 0.0;
  std::vector<int> seen;
  for (int b = 0; b < nb; ++b) if (std::find(seen.begin(), seen.end(), bstr[b]) == seen.end()) seen.push_back(bstr[b]);
  for (int s : seen) {
    double prev = 0.0;
    for (int b = 0; b < nb; ++b) {
      if (bstr[b] != s) continue;
      const double inc = bcum[b] - prev; prev = bcum[b];
      const double* hit = std::lower_bound(tj, tj + m, btime[b]);
      if (hit == tj + m || *hit != btime[b]) return -1;
      const std::size_t j = static_cast<std::size_t>(hit - tj);
      for (int i = 0; i < nr; ++i) {
        if (rstr[i] != s) continue;
        const double lp = lp_by_time ? LP[static_cast<std::size_t>(i) * m + j] : LP[i];
        out[static_cast<std::size_t>(i) * m + j] = inc * std::exp(lp);
      }
    }
  }
  return 0;
} catch (...) { return -9; }


// covariates that are linear combinations of the others and of the constant (survcore/scaling.hpp):
// out (room for p) receives their 0-based indices, *ndep their number
int sc_collinear(int n, int p, const double* X, double tol, int* out, int* ndep) try {
  if (n < 0 || p < 0 || !ndep || (n > 0 && p > 0 && (!X || !out)) || !(tol > 0)) return -1;
  for (std::size_t k = 0; k < static_cast<std::size_t>(n) * p; ++k) if (!std::isfinite(X[k])) return -2;
  const std::vector<int> d = collinear_columns(static_cast<std::size_t>(n), static_cast<std::size_t>(p), X, tol);
  *ndep = static_cast<int>(d.size());
  for (std::size_t k = 0; k < d.size(); ++k) out[k] = d[k];
  return 0;
} catch (...) { return -9; }

}  // extern "C"
