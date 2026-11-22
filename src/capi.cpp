#include <type_traits>
#include <cstdint>
#include "survcore/survcore.h"
#include <cmath>
#include <limits>
#include <vector>
#include "survcore/compare.hpp"
#include "survcore/cox.hpp"
#include "survcore/directcif.hpp"
#include "survcore/directcif_ext.hpp"
#include "survcore/metrics.hpp"
#include "survcore/finegray.hpp"
#include "survcore/nonparametric.hpp"
#include "survcore/tests.hpp"
#include "survcore/scaling.hpp"

using namespace survcore;

namespace {
bool finite_all(const double* x, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
  return true;
}
// finite and not negative (times, weights)
bool nonneg_all(const double* x, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) if (!std::isfinite(x[i]) || x[i] < 0.0) return false;
  return true;
}
int check_common(int n, int p, const double* time, const double* X) {
  if (n <= 0 || p < 0 || !time || (p > 0 && !X)) return -1;
  if (!finite_all(time, n) || !finite_all(X, static_cast<std::size_t>(n) * p)) return -2;
  for (int i = 0; i < n; ++i) if (time[i] < 0) return -2;
  return 0;
}
void export_fit(const FitResult& r, std::size_t p, int n_events, double* beta, double* var, sc_fit_info* info) {
  for (std::size_t j = 0; j < p; ++j) beta[j] = r.beta[j];
  if (var)
    for (std::size_t j = 0; j < p * p; ++j)
      var[j] = r.var.n == p ? r.var.a[j] : std::numeric_limits<double>::quiet_NaN();
  if (info) {
    info->loglik_init = r.ll_init; info->loglik = r.ll; info->grad_maxabs = r.grad_maxabs;
    info->iterations = r.iterations; info->converged = r.converged; info->flags = r.flags;
    info->n_events = n_events;
  }
}
// Monotone likelihood (Heinze & Schemper 2001, Biometrics 57:114-119): the
// log-likelihood converges while a coefficient keeps moving towards infinity.
// Flag when the last Newton step still moved the linear predictor by more than
// 0.05 per covariate standard deviation, or when the fitted effect is absurdly large.
void flag_monotone(FitResult& r, int n, int p, const double* X) {
  for (int j = 0; j < p; ++j) {
    double m = 0, s = 0;
    for (int i = 0; i < n; ++i) m += X[static_cast<std::size_t>(i) * p + j];
    m /= n;
    for (int i = 0; i < n; ++i) { const double d = X[static_cast<std::size_t>(i) * p + j] - m; s += d * d; }
    s = std::sqrt(s / n);
    const double last = r.last_step.empty() ? 0.0 : std::fabs(r.last_step[j]);
    if (last * s > 0.05 || std::fabs(r.beta[j]) * 2.0 * s > 30.0) r.flags |= FIT_POSSIBLE_MONOTONE;
  }
}
void export_eval(const Eval& ev, std::size_t p, double* ll, double* g, double* H) {
  if (ll) *ll = ev.ll;
  if (g) for (std::size_t j = 0; j < p; ++j) g[j] = ev.g[j];
  if (H) for (std::size_t j = 0; j < p * p; ++j) H[j] = ev.H.a[j];
}
}  // namespace

extern "C" {

int sc_kaplan_meier(int n, const double* time, const int* status, const double* weights,
                    double* out_time, double* out_nrisk, double* out_nevent,
                    double* out_surv, double* out_se, int* m) try {
  if (n <= 0 || !time || !status || !out_time || !out_surv || !m) return -1;
  if (!nonneg_all(time, n) || (weights && !nonneg_all(weights, n))) return -2;
  std::vector<double> t(time, time + n), w;
  if (weights) w.assign(weights, weights + n);
  std::vector<int> s(status, status + n);
  KMResult r = kaplan_meier(t, s, w);
  *m = static_cast<int>(r.time.size());
  for (std::size_t i = 0; i < r.time.size(); ++i) {
    out_time[i] = r.time[i]; out_surv[i] = r.surv[i];
    if (out_nrisk) out_nrisk[i] = r.n_risk[i];
    if (out_nevent) out_nevent[i] = r.n_event[i];
    if (out_se) out_se[i] = r.std_err[i];
  }
  return 0;
} catch (...) { return -9; }

int sc_nelson_aalen(int n, const double* time, const int* status, const double* weights,
                    double* out_time, double* out_cumhaz, double* out_se, int* m) try {
  if (n <= 0 || !time || !status || !out_time || !out_cumhaz || !m) return -1;
  if (!nonneg_all(time, n) || (weights && !nonneg_all(weights, n))) return -2;
  std::vector<double> w;
  if (weights) w.assign(weights, weights + n);
  KMResult r = kaplan_meier(std::vector<double>(time, time + n), std::vector<int>(status, status + n), w);
  *m = static_cast<int>(r.time.size());
  for (std::size_t i = 0; i < r.time.size(); ++i) { out_time[i] = r.time[i]; out_cumhaz[i] = r.cumhaz[i]; if (out_se) out_se[i] = r.cumhaz_se[i]; }
  return 0;
} catch (...) { return -9; }

int sc_rmst(int n, const double* time, const int* status, double tau, double* mean, double* se) try {
  if (n <= 0 || !time || !status || !mean || !se) return -1;
  if (!nonneg_all(time, n) || !(tau > 0.0) || !std::isfinite(tau)) return -2;
  KMResult r = kaplan_meier(std::vector<double>(time, time + n), std::vector<int>(status, status + n), {});
  rmst(r, tau, *mean, *se);
  return 0;
} catch (...) { return -9; }

int sc_aalen_johansen(int n, const double* time, const int* cause, int K,
                      double* out_time, double* out_surv, double* out_cif, double* out_se, int* m) try {
  if (n <= 0 || K <= 0 || !time || !cause || !out_time || !out_cif || !m) return -1;
  if (!nonneg_all(time, n)) return -2;
  std::vector<double> t(time, time + n);
  std::vector<int> c(cause, cause + n);
  AJResult r = aalen_johansen(t, c, K);
  *m = static_cast<int>(r.time.size());
  for (std::size_t i = 0; i < r.time.size(); ++i) {
    out_time[i] = r.time[i];
    if (out_surv) out_surv[i] = r.surv[i];
  }
  for (std::size_t i = 0; i < r.cif.size(); ++i) { out_cif[i] = r.cif[i]; if (out_se) out_se[i] = r.std_err[i]; }
  return 0;
} catch (...) { return -9; }

namespace {
int check_ext(int n, const double* entry, const double* time, const int* strata, const double* weights) {
  if (entry) {
    if (!finite_all(entry, n)) return -2;
    for (int i = 0; i < n; ++i) if (!(entry[i] < time[i])) return -6;
  }
  if (strata) for (int i = 0; i < n; ++i) if (strata[i] < 0) return -1;
  if (weights) { if (!finite_all(weights, n)) return -2; for (int i = 0; i < n; ++i) if (weights[i] < 0) return -2; }
  return 0;
}
}  // namespace

int sc_coxph_ext(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                 const double* X, const double* weights, int ties, int max_iter, double tol,
                 double* beta, double* var, sc_fit_info* info) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!status || !beta || p == 0) return -1;
  if (int rc = check_ext(n, entry, time, strata, weights)) return rc;
  const ColumnScale cs = scale_columns(static_cast<std::size_t>(n), static_cast<std::size_t>(p), X);   // scaling.hpp
  const double* Xs = cs.data(X);
  CoxModel m(n, p, time, status, Xs, weights, ties == 1 ? EFRON : BRESLOW, entry, strata);
  if (m.n_events() == 0) return -3;
  std::vector<double> b0(beta, beta + p);
  coef_to_scaled(cs, 0, b0.data());
  FitResult r = newton_fit(m, b0, max_iter, tol);
  flag_monotone(r, n, p, Xs);
  export_fit(r, p, static_cast<int>(m.n_events()), beta, var, info);
  coef_from_scaled(cs, 0, beta);
  var_from_scaled(cs, 0, static_cast<std::size_t>(p), var);
  return 0;
} catch (...) { return -9; }

int sc_coxph(int n, int p, const double* time, const int* status, const double* X,
             const double* weights, int ties, int max_iter, double tol,
             double* beta, double* var, sc_fit_info* info) {
  return sc_coxph_ext(n, p, nullptr, time, status, nullptr, X, weights, ties, max_iter, tol, beta, var, info);
}

int sc_coxph_ext_baseline(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                          const double* X, const double* weights, int ties, const double* beta, const double* zref,
                          double* out_time, double* out_cumhaz, int* out_strata, int* m) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!status || !beta || !out_time || !out_cumhaz || !m) return -1;
  if (int rc = check_ext(n, entry, time, strata, weights)) return rc;
  CoxModel mod(n, p, time, status, X, weights, ties == 1 ? EFRON : BRESLOW, entry, strata);
  std::vector<double> t, H;
  std::vector<int> s;
  mod.baseline(std::vector<double>(beta, beta + p), zref, t, H, &s);
  *m = static_cast<int>(t.size());
  for (std::size_t i = 0; i < t.size(); ++i) {
    out_time[i] = t[i]; out_cumhaz[i] = H[i];
    if (out_strata) out_strata[i] = s[i];
  }
  return 0;
} catch (...) { return -9; }

int sc_coxph_resid(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                   const double* X, const double* weights, int ties, const double* beta, double* mart, double* xbar) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!status || !beta || !mart || !xbar || p == 0) return -1;
  if (int rc = check_ext(n, entry, time, strata, weights)) return rc;
  CoxModel mod(n, p, time, status, X, weights, ties == 1 ? EFRON : BRESLOW, entry, strata);
  std::vector<double> u, m, xb;
  mod.score_residuals(std::vector<double>(beta, beta + p), u, &m, &xb);
  for (int i = 0; i < n; ++i) mart[i] = m[i];
  for (std::size_t i = 0; i < xb.size(); ++i) xbar[i] = xb[i];
  return 0;
} catch (...) { return -9; }

int sc_coxph_score(int n, int p, const double* entry, const double* time, const int* status, const int* strata,
                   const double* X, const double* weights, int ties, const double* beta, double* U) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!status || !beta || !U || p == 0) return -1;
  if (int rc = check_ext(n, entry, time, strata, weights)) return rc;
  CoxModel mod(n, p, time, status, X, weights, ties == 1 ? EFRON : BRESLOW, entry, strata);
  std::vector<double> u;
  mod.score_residuals(std::vector<double>(beta, beta + p), u);
  for (std::size_t i = 0; i < u.size(); ++i) U[i] = u[i];
  return 0;
} catch (...) { return -9; }

int sc_coxph_eval(int n, int p, const double* time, const int* status, const double* X,
                  const double* weights, int ties, const double* beta,
                  double* loglik, double* grad, double* info_matrix) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!status || !beta) return -1;
  CoxModel m(n, p, time, status, X, weights, ties == 1 ? EFRON : BRESLOW);
  Eval ev;
  m.eval(std::vector<double>(beta, beta + p), ev, true);
  export_eval(ev, p, loglik, grad, info_matrix);
  return 0;
} catch (...) { return -9; }

int sc_finegray(int n, int p, const double* time, const int* cause, const double* X,
                int interest, int max_iter, double tol,
                double* beta, double* var_model, double* var_robust, sc_fit_info* info) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause || !beta || p == 0 || interest == 0) return -1;
  const ColumnScale cs = scale_columns(static_cast<std::size_t>(n), static_cast<std::size_t>(p), X);   // scaling.hpp
  const double* Xs = cs.data(X);
  FineGrayModel m(n, p, time, cause, Xs, interest);
  if (m.n_events() == 0) return -3;
  std::vector<double> b0(beta, beta + p);
  coef_to_scaled(cs, 0, b0.data());
  FitResult r = newton_fit(m, b0, max_iter, tol);
  flag_monotone(r, n, p, Xs);
  export_fit(r, p, static_cast<int>(m.n_events()), beta, var_model, info);
  coef_from_scaled(cs, 0, beta);
  var_from_scaled(cs, 0, static_cast<std::size_t>(p), var_model);
  if (var_robust) {
    const std::size_t P = p;
    if (r.var.n == P) {
      Mat meat;
      m.robust_meat(r.beta, meat);
      for (std::size_t a = 0; a < P; ++a)
        for (std::size_t b = 0; b < P; ++b) {
          double s = 0.0;
          for (std::size_t u = 0; u < P; ++u)
            for (std::size_t v = 0; v < P; ++v) s += r.var(a, u) * meat(u, v) * r.var(v, b);
          var_robust[a * P + b] = s;
        }
    } else {
      for (std::size_t j = 0; j < P * P; ++j) var_robust[j] = std::numeric_limits<double>::quiet_NaN();
    }
    var_from_scaled(cs, 0, P, var_robust);
  }
  return 0;
} catch (...) { return -9; }

int sc_finegray_eval(int n, int p, const double* time, const int* cause, const double* X,
                     int interest, const double* beta,
                     double* loglik, double* grad, double* info_matrix) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause || !beta || interest == 0) return -1;
  FineGrayModel m(n, p, time, cause, X, interest);
  Eval ev;
  m.eval(std::vector<double>(beta, beta + p), ev, true);
  export_eval(ev, p, loglik, grad, info_matrix);
  return 0;
} catch (...) { return -9; }

int sc_coxph_baseline(int n, int p, const double* time, const int* status, const double* X,
                      const double* weights, int ties, const double* beta, const double* zref,
                      double* out_time, double* out_cumhaz, int* m) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!status || !beta || !out_time || !out_cumhaz || !m) return -1;
  CoxModel mod(n, p, time, status, X, weights, ties == 1 ? EFRON : BRESLOW);
  std::vector<double> t, H;
  mod.baseline(std::vector<double>(beta, beta + p), zref, t, H);
  *m = static_cast<int>(t.size());
  for (std::size_t i = 0; i < t.size(); ++i) { out_time[i] = t[i]; out_cumhaz[i] = H[i]; }
  return 0;
} catch (...) { return -9; }

int sc_finegray_baseline(int n, int p, const double* time, const int* cause, const double* X,
                         int interest, const double* beta, const double* zref,
                         double* out_time, double* out_cumhaz, int* m) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause || !beta || !out_time || !out_cumhaz || !m || interest == 0) return -1;
  FineGrayModel mod(n, p, time, cause, X, interest);
  std::vector<double> H;
  mod.baseline(std::vector<double>(beta, beta + p), zref, H);
  *m = static_cast<int>(H.size());
  for (std::size_t i = 0; i < H.size(); ++i) { out_time[i] = mod.event_times()[i]; out_cumhaz[i] = H[i]; }
  return 0;
} catch (...) { return -9; }

int sc_predict_from_cumhaz(int m, const double* cumhaz, int n_new, const double* lp, int mode, double* out) try {
  if (m < 0 || n_new < 0 || !cumhaz || !lp || !out) return -1;
  for (int i = 0; i < n_new; ++i) {
    const double r = std::exp(lp[i]);
    for (int j = 0; j < m; ++j) {
      const double a = -cumhaz[j] * r;
      out[static_cast<std::size_t>(i) * m + j] = mode == 1 ? -std::expm1(a) : std::exp(a);
    }
  }
  return 0;
} catch (...) { return -9; }

int sc_csc_cif(int n, int p, const double* time, const int* cause, const double* X, int K,
               const double* betas, int ties, int n_new, const double* Xnew, int product_limit,
               double* out_time, double* out_cif, int* m) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause || !betas || K <= 0 || n_new < 0 || (n_new > 0 && !Xnew) || !out_time || !out_cif || !m) return -1;
  const std::size_t P = p;
  std::vector<double> zref(P, 0.0);                 // centre at the column means for stability
  for (int i = 0; i < n; ++i) for (std::size_t j = 0; j < P; ++j) zref[j] += X[static_cast<std::size_t>(i) * P + j] / n;
  std::vector<std::vector<double>> tk(K), Hk(K);
  std::vector<int> st(n);
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < n; ++i) st[i] = cause[i] == k + 1;
    CoxModel mod(n, p, time, st.data(), X, nullptr, ties == 1 ? EFRON : BRESLOW);
    mod.baseline(std::vector<double>(betas + k * P, betas + (k + 1) * P), zref.data(), tk[k], Hk[k]);
  }
  std::vector<double> grid;
  for (int k = 0; k < K; ++k) grid.insert(grid.end(), tk[k].begin(), tk[k].end());
  std::sort(grid.begin(), grid.end());
  grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
  const std::size_t M = grid.size();
  *m = static_cast<int>(M);
  std::vector<double> dH(M * K, 0.0);               // hazard increments on the common grid
  for (int k = 0; k < K; ++k) {
    std::size_t q = 0; double prev = 0.0;
    for (std::size_t j = 0; j < M && q < tk[k].size(); ++j)
      if (tk[k][q] == grid[j]) { dH[j * K + k] = Hk[k][q] - prev; prev = Hk[k][q]; ++q; }
  }
  for (std::size_t j = 0; j < M; ++j) out_time[j] = grid[j];
  std::vector<double> rr(K), F(K);
  for (int i = 0; i < n_new; ++i) {
    for (int k = 0; k < K; ++k) {
      double lp = 0.0;
      for (std::size_t j = 0; j < P; ++j) lp += (Xnew[static_cast<std::size_t>(i) * P + j] - zref[j]) * betas[k * P + j];
      rr[k] = std::exp(lp);
    }
    double S = 1.0;
    std::fill(F.begin(), F.end(), 0.0);
    for (std::size_t j = 0; j < M; ++j) {
      double tot = 0.0;
      for (int k = 0; k < K; ++k) { const double hz = dH[j * K + k] * rr[k]; F[k] += S * hz; tot += hz; }
      S *= product_limit ? std::fmax(0.0, 1.0 - tot) : std::exp(-tot);
      for (int k = 0; k < K; ++k) out_cif[(static_cast<std::size_t>(i) * M + j) * K + k] = F[k];
    }
  }
  return 0;
} catch (...) { return -9; }

namespace {
int check_groups(int n, const double* time, const int* ev, const int* group, int K) {
  if (n <= 0 || K < 2 || !time || !ev || !group) return -1;
  if (!finite_all(time, n)) return -2;
  for (int i = 0; i < n; ++i) if (group[i] < 0 || group[i] >= K) return -2;
  return 0;
}
void export_test(const KSampleResult& r, int K, double* stat, int* df, double* p, double* score, double* var) {
  if (stat) *stat = r.statistic;
  if (df) *df = r.df;
  if (p) *p = r.p_value;
  if (score) for (int g = 0; g < K; ++g) score[g] = r.score[g];
  if (var) for (int g = 0; g < K * K; ++g) var[g] = r.var.a[g];
}
}  // namespace

int sc_logrank(int n, const double* time, const int* status, const int* group, int K,
               const int* strata, double rho, double* stat, int* df, double* p_value,
               double* out_obs, double* out_exp, double* out_var) try {
  if (int rc = check_groups(n, time, status, group, K)) return rc;
  KSampleResult r = logrank_test(std::vector<double>(time, time + n), std::vector<int>(status, status + n),
                                 std::vector<int>(group, group + n), K, strata, rho);
  export_test(r, K, stat, df, p_value, nullptr, out_var);
  for (int g = 0; g < K; ++g) { if (out_obs) out_obs[g] = r.observed[g]; if (out_exp) out_exp[g] = r.expected[g]; }
  return 0;
} catch (...) { return -9; }

int sc_gray_test(int n, const double* time, const int* cause, const int* group, int K,
                 const int* strata, int interest, double rho, double* stat, int* df,
                 double* p_value, double* out_score, double* out_var) try {
  if (int rc = check_groups(n, time, cause, group, K)) return rc;
  if (interest == 0) return -1;
  KSampleResult r = gray_test(std::vector<double>(time, time + n), std::vector<int>(cause, cause + n),
                              std::vector<int>(group, group + n), K, strata, interest, rho);
  export_test(r, K, stat, df, p_value, out_score, out_var);
  return 0;
} catch (...) { return -9; }

double sc_chisq_upper(double x, double df) try { return chisq_upper(x, df); }
catch (...) { return std::numeric_limits<double>::quiet_NaN(); }

}  // extern "C"

// ---------------- Phase 2 ----------------
namespace {
int make_spec(int p, int K, const int* alpha_free, const double* alpha_value, const int* baseline,
              const int* n_knots, const double* knots, DirectSpec& s) {
  if (p < 0 || K < 1 || !alpha_free) return -1;
  s.K = K; s.p = static_cast<std::size_t>(p);
  s.alpha_free.assign(alpha_free, alpha_free + K);
  s.alpha_value.assign(K, 0.0);
  if (alpha_value) s.alpha_value.assign(alpha_value, alpha_value + K);
  s.baseline.assign(K, BASE_GOMPERTZ);
  if (baseline) s.baseline.assign(baseline, baseline + K);
  s.knots.assign(K, std::vector<double>());
  std::size_t pos = 0;
  for (int k = 0; k < K; ++k) {
    if (!std::isfinite(s.alpha_value[k])) return -2;
    if (s.baseline[k] < 0 || s.baseline[k] >= BASE_COUNT) return -1;
    const int nk = n_knots ? n_knots[k] : 0;
    if (nk < 0 || (nk > 0 && !knots)) return -1;
    if (s.baseline[k] == BASE_SPLINE && nk > 0) {
      if (nk < 3) return -1;                                   // boundary knots plus at least one interior knot
      s.knots[k].assign(knots + pos, knots + pos + nk);
      for (int j = 0; j < nk; ++j) if (!std::isfinite(s.knots[k][j]) || (j > 0 && !(s.knots[k][j] > s.knots[k][j - 1]))) return -2;
    }
    pos += static_cast<std::size_t>(nk);
  }
  return 0;
}
int check_direct(int n, int p, int K, const double* time, const int* cause, const double* X, const DirectSpec* s = nullptr) {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause) return -1;
  if (s) for (int k = 0; k < K; ++k) if (base_needs_positive_time(s->baseline[k])) for (int i = 0; i < n; ++i) if (!(time[i] > 0.0)) return -2;
  std::vector<int> cnt(K + 1, 0);
  for (int i = 0; i < n; ++i) { if (cause[i] < 0 || cause[i] > K) return -2; ++cnt[cause[i]]; }
  for (int k = 1; k <= K; ++k) if (cnt[k] == 0) return -3;
  return 0;
}
int check_pred(int n, const double* time, const int* cause, const double* pred) {
  if (n <= 0 || !time || !cause || !pred) return -1;
  if (!finite_all(time, n) || !finite_all(pred, n)) return -2;
  return 0;
}
}  // namespace

extern "C" {

int sc_directcif_nparam(int p, int K, const int* alpha_free, const int* baseline, const int* n_knots) try {
  if (K < 1) return -1;
  DirectSpec s;
  std::vector<double> dummy;
  if (n_knots) { int tot = 0; for (int k = 0; k < K; ++k) tot += n_knots[k] > 0 ? n_knots[k] : 0; dummy.resize(tot); for (int j = 0; j < tot; ++j) dummy[j] = j; }
  if (make_spec(p, K, alpha_free, nullptr, baseline, n_knots, dummy.empty() ? nullptr : dummy.data(), s)) return -1;
  return static_cast<int>(s.nparam());
} catch (...) { return -9; }

int sc_directcif_default_knots(int n, const double* time, const int* cause, int k, int df, double* out_knots, int* n_knots) try {
  if (n <= 0 || !time || !cause || !out_knots || !n_knots || k < 1 || df < 1) return -1;
  const std::vector<double> kn = direct_default_knots(n, time, cause, k, df);
  *n_knots = static_cast<int>(kn.size());
  for (std::size_t j = 0; j < kn.size(); ++j) out_knots[j] = kn[j];
  return 0;
} catch (...) { return -9; }

const char* sc_directcif_base_param_name(int baseline, int j) {
  static const char* nm[BASE_COUNT][3] = {{"rho", "log_tau", ""}, {"rho", "log_tau", "zeta"}, {"log_lambda", "log_kappa", ""},
    {"logit_p", "log_lambda", "log_kappa"}, {"logit_p", "log_lambda", "log_kappa"}, {"logit_p", "mu", "log_sigma"},
    {"log_b", "c", "logit_p"}, {"gamma", "gamma", "gamma"}};
  if (baseline < 0 || baseline >= BASE_COUNT || j < 0) return "";
  return baseline == BASE_SPLINE ? "gamma" : (j < 3 ? nm[baseline][j] : "");
}

int sc_directcif_fit(int n, int p, int K, const double* time, const int* cause, const double* X,
                     const int* alpha_free, const double* alpha_value, const int* baseline,
                     const int* n_knots, const double* knots, int constraint,
                     const double* init, int max_iter, double tol,
                     double* theta, double* var, sc_fit_info* info, double* diag) try {
  DirectSpec s;
  if (int rc = make_spec(p, K, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (int rc = check_direct(n, p, K, time, cause, X, &s)) return rc;
  if (!theta) return -1;
  const std::size_t P = s.nparam();
  DirectCIFModel m(n, time, cause, X, s);
  std::vector<double> st;
  if (init) { if (!finite_all(init, P)) return -2; st.assign(init, init + P); }
  const DirectFit f = direct_fit(m, n, time, cause, init ? &st : nullptr, constraint != 0, max_iter, tol);
  int ne = 0;
  for (int i = 0; i < n; ++i) ne += cause[i] != 0;
  export_fit(f.fit, P, ne, theta, var, info);
  if (diag) {
    diag[0] = f.diag.n_exceed_range; diag[1] = f.diag.max_sum_range; diag[2] = f.diag.n_exceed_inf;
    diag[3] = f.diag.max_sum_inf; diag[4] = f.diag.first_crossing_time;
  }
  return 0;
} catch (...) { return -9; }

int sc_directcif_fit_ic(int n, int p, int K, const double* time, const double* time2, const int* cause, const double* X,
                        const int* alpha_free, const double* alpha_value, const int* baseline,
                        const int* n_knots, const double* knots, int constraint,
                        const double* init, int max_iter, double tol,
                        double* theta, double* var, sc_fit_info* info, double* diag) try {
  DirectSpec s;
  if (int rc = make_spec(p, K, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (!time || !time2 || !cause) return -1;
  // midpoints of the intervals: the data-driven start and the checks of the common entry point
  std::vector<double> mid(n);
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(time[i]) || time[i] < 0.0) return -2;
    if (cause[i] > 0) {
      if (!std::isfinite(time2[i]) || time2[i] < time[i] || !(time2[i] > 0.0)) return -2;
      mid[i] = 0.5 * (time[i] + time2[i]);
    } else {
      if (!(time[i] > 0.0)) return -2;
      mid[i] = time[i];
    }
  }
  if (int rc = check_direct(n, p, K, mid.data(), cause, X, &s)) return rc;
  if (!theta) return -1;
  const std::size_t P = s.nparam();
  DirectCIFModel m(n, time, cause, X, s);
  m.set_time2(time2);
  std::vector<double> st;
  if (init) { if (!finite_all(init, P)) return -2; st.assign(init, init + P); }
  const DirectFit f = direct_fit(m, n, mid.data(), cause, init ? &st : nullptr, constraint != 0, max_iter, tol);
  int ne = 0;
  for (int i = 0; i < n; ++i) ne += cause[i] != 0;
  export_fit(f.fit, P, ne, theta, var, info);
  if (diag) {
    diag[0] = f.diag.n_exceed_range; diag[1] = f.diag.max_sum_range; diag[2] = f.diag.n_exceed_inf;
    diag[3] = f.diag.max_sum_inf; diag[4] = f.diag.first_crossing_time;
  }
  return 0;
} catch (...) { return -9; }

int sc_directcif_eval_ic(int n, int p, int K, const double* time, const double* time2, const int* cause, const double* X,
                         const int* alpha_free, const double* alpha_value, const int* baseline,
                         const int* n_knots, const double* knots, const double* theta,
                         double* loglik, double* grad, double* info_matrix) try {
  DirectSpec s;
  if (int rc = make_spec(p, K, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (!time || !time2 || !cause || !theta) return -1;
  DirectCIFModel m(n, time, cause, X, s);
  m.set_time2(time2);
  Eval ev;
  m.eval(std::vector<double>(theta, theta + s.nparam()), ev, true);
  export_eval(ev, s.nparam(), loglik, grad, info_matrix);
  return 0;
} catch (...) { return -9; }

int sc_directcif_eval(int n, int p, int K, const double* time, const int* cause, const double* X,
                      const int* alpha_free, const double* alpha_value, const int* baseline,
                      const int* n_knots, const double* knots, const double* theta,
                      double* loglik, double* grad, double* info_matrix) try {
  DirectSpec s;
  if (int rc = make_spec(p, K, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (int rc = check_direct(n, p, K, time, cause, X, &s); rc != 0 && rc != -3) return rc;
  if (!cause || !theta) return -1;
  DirectCIFModel m(n, time, cause, X, s);
  Eval ev;
  m.eval(std::vector<double>(theta, theta + s.nparam()), ev, true);
  export_eval(ev, s.nparam(), loglik, grad, info_matrix);
  return 0;
} catch (...) { return -9; }

int sc_directcif_predict(int p, int K, const int* alpha_free, const double* alpha_value, const int* baseline,
                         const int* n_knots, const double* knots,
                         const double* theta, const double* var, int k, int n_new, const double* Xnew,
                         int nt, const double* times, double* out_cif, double* out_se) try {
  DirectSpec s;
  if (int rc = make_spec(p, K, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (!theta || !times || !out_cif || n_new <= 0 || nt <= 0 || k < 1 || k > K || (p > 0 && !Xnew)) return -1;
  if (out_se && !var) return -1;
  const std::size_t P = s.nparam();
  std::vector<double> dummy_t(1, 0.0); std::vector<int> dummy_c(1, 0);
  DirectCIFModel m(0, dummy_t.data(), dummy_c.data(), Xnew, s);
  const std::vector<double> th(theta, theta + P);
  std::vector<double> g;
  for (int i = 0; i < n_new; ++i)
    for (int t = 0; t < nt; ++t) {
      const double F = m.cif_grad(th, k - 1, Xnew + static_cast<std::size_t>(i) * p, times[t], g);
      out_cif[static_cast<std::size_t>(i) * nt + t] = F;
      if (out_se) {
        double v = 0.0;
        for (std::size_t a = 0; a < P; ++a) { if (g[a] == 0.0) continue; for (std::size_t b = 0; b < P; ++b) v += g[a] * var[a * P + b] * g[b]; }
        out_se[static_cast<std::size_t>(i) * nt + t] = v >= 0.0 ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN();
      }
    }
  return 0;
} catch (...) { return -9; }

int sc_brier(int n, const double* time, const int* cause, int interest, const double* pred, double t,
             double* brier, double* brier_null, double* ipa) try {
  if (int rc = check_pred(n, time, cause, pred)) return rc;
  const BrierResult r = brier_score(std::vector<double>(time, time + n), std::vector<int>(cause, cause + n), interest,
                                    std::vector<double>(pred, pred + n), t);
  if (brier) *brier = r.brier;
  if (brier_null) *brier_null = r.brier_null;
  if (ipa) *ipa = r.ipa;
  return 0;
} catch (...) { return -9; }

int sc_auc(int n, const double* time, const int* cause, int interest, const double* pred, double t, double* auc) try {
  if (int rc = check_pred(n, time, cause, pred)) return rc;
  if (!auc) return -1;
  *auc = auc_t(std::vector<double>(time, time + n), std::vector<int>(cause, cause + n), interest,
               std::vector<double>(pred, pred + n), t);
  return 0;
} catch (...) { return -9; }

int sc_calibration(int n, const double* time, const int* cause, int interest, const double* pred, double t,
                   int groups, double* out_mean_pred, double* out_observed, double* out_size, double* error, int* m) try {
  if (int rc = check_pred(n, time, cause, pred)) return rc;
  if (groups < 1 || !out_mean_pred || !out_observed || !m) return -1;
  const CalibrationResult r = calibration(std::vector<double>(time, time + n), std::vector<int>(cause, cause + n),
                                          interest, std::vector<double>(pred, pred + n), t, groups);
  *m = static_cast<int>(r.mean_pred.size());
  for (std::size_t g = 0; g < r.mean_pred.size(); ++g) {
    out_mean_pred[g] = r.mean_pred[g]; out_observed[g] = r.observed[g];
    if (out_size) out_size[g] = r.size[g];
  }
  if (error) *error = r.error;
  return 0;
} catch (...) { return -9; }

int sc_compare_models2(int n, int p, int K, const double* time, const int* cause, const double* X,
                       int interest, int nt, const double* times, int folds, unsigned long long seed,
                       int n_direct, const int* baseline, const int* link, const int* spline_df,
                       const int* estimation, const int* cond, double* out_metrics, double* out_fit) try {
  if (int rc = check_direct(n, p, K, time, cause, X)) return rc;
  if (p < 1 || nt < 1 || !times || !out_metrics || interest < 1 || interest > K || folds > n || n_direct < 0) return -1;
  if (n_direct > 0 && (!baseline || !link)) return -1;
  std::vector<DirectChoice> ch(n_direct);
  for (int j = 0; j < n_direct; ++j) {
    if (baseline[j] < 0 || baseline[j] >= BASE_COUNT || link[j] < 0 || link[j] > 2) return -1;
    ch[j].baseline = baseline[j]; ch[j].link = link[j]; ch[j].spline_df = spline_df ? spline_df[j] : 3;
    ch[j].estimation = estimation ? estimation[j] : 0; ch[j].cond = cond ? cond[j] : COND_LOGISTIC;
    if (ch[j].estimation < 0 || ch[j].estimation > 2 || ch[j].cond < 0 || ch[j].cond >= COND_COUNT) return -1;
    if (ch[j].estimation == 2 && (K != 2 || baseline[j] == BASE_WEIBULL || baseline[j] == BASE_SPLINE)) return -1;
    if (base_needs_positive_time(baseline[j]) || (ch[j].estimation == 2 && ch[j].cond != COND_LOGISTIC))
      for (int i = 0; i < n; ++i) if (!(time[i] > 0.0)) return -2;
  }
  const CompareResult r = compare_models(ch, n, p, K, time, cause, X, interest, std::vector<double>(times, times + nt), folds, seed);
  for (std::size_t i = 0; i < r.metrics.size(); ++i) out_metrics[i] = r.metrics[i];
  if (out_fit) for (std::size_t i = 0; i < r.fit.size(); ++i) out_fit[i] = r.fit[i];
  return 0;
} catch (...) { return -9; }

int sc_compare_models(int n, int p, int K, const double* time, const int* cause, const double* X,
                      int interest, int nt, const double* times, int folds, unsigned long long seed,
                      int n_direct, const int* baseline, const int* link, const int* spline_df,
                      double* out_metrics, double* out_fit) {
  return sc_compare_models2(n, p, K, time, cause, X, interest, nt, times, folds, seed, n_direct, baseline, link, spline_df,
                            nullptr, nullptr, out_metrics, out_fit);
}

int sc_cv_folds(int n, int folds, unsigned long long seed, int* out_fold) try {
  if (n <= 0 || folds < 1 || !out_fold) return -1;
  const std::vector<int> f = cv_folds(n, folds, seed);
  for (int i = 0; i < n; ++i) out_fold[i] = f[i];
  return 0;
} catch (...) { return -9; }

}  // extern "C"

namespace {
int one_spec(int p, int alpha_free, double alpha_value, int baseline, int n_knots, const double* knots, DirectSpec& s) {
  return make_spec(p, 1, &alpha_free, &alpha_value, &baseline, &n_knots, knots, s);
}
int check_subdist(int n, int p, const double* time, const int* cause, const double* X, int interest, const DirectSpec& s) {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause || interest < 1) return -1;
  int ne = 0;
  for (int i = 0; i < n; ++i) { if (cause[i] < 0) return -2; ne += cause[i] == interest; }
  if (base_needs_positive_time(s.baseline[0])) for (int i = 0; i < n; ++i) if (!(time[i] > 0.0)) return -2;
  return ne == 0 ? -3 : 0;
}
int check_constrained(int n, int p, const double* time, const int* cause, const double* X, int family, const DirectSpec& s) {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause || family < 0 || family >= COND_COUNT) return -1;
  if (s.baseline[0] == BASE_WEIBULL || s.baseline[0] == BASE_SPLINE) return -1;   // cause 1 needs a plateau
  int n1 = 0, n2 = 0;
  for (int i = 0; i < n; ++i) { if (cause[i] < 0 || cause[i] > 2) return -2; n1 += cause[i] == 1; n2 += cause[i] == 2; }
  if (family != COND_LOGISTIC || base_needs_positive_time(s.baseline[0])) for (int i = 0; i < n; ++i) if (!(time[i] > 0.0)) return -2;
  return (n1 == 0 || n2 == 0) ? -3 : 0;
}
int check_constrainedK(int n, int p, const double* time, const int* cause, const double* X, int K, const int* family,
                       const DirectSpec& s) {
  if (int rc = check_common(n, p, time, X)) return rc;
  if (!cause || !family || K < 2) return -1;
  if (s.baseline[0] == BASE_WEIBULL || s.baseline[0] == BASE_SPLINE) return -1;
  bool need_pos = base_needs_positive_time(s.baseline[0]);
  for (int k = 0; k < K - 1; ++k) { if (family[k] < 0 || family[k] >= COND_COUNT) return -1; if (family[k] != COND_LOGISTIC) need_pos = true; }
  std::vector<int> cnt(static_cast<std::size_t>(K + 1), 0);
  for (int i = 0; i < n; ++i) { if (cause[i] < 0 || cause[i] > K) return -2; cnt[static_cast<std::size_t>(cause[i])] += 1; }
  if (need_pos) for (int i = 0; i < n; ++i) if (!(time[i] > 0.0)) return -2;
  for (int k = 1; k <= K; ++k) if (cnt[static_cast<std::size_t>(k)] == 0) return -3;
  return 0;
}
// Covariates of the time distributions G_k (sc_constrainedK_fit_cov). wuse, vuse: (K - 1) x nw and
// (K - 1) x nv flags (row k - 2 for cause k) of the columns of W and V used for cause k.
int make_cond_cov(int n_rows, int K, int nw, const double* W, const int* wuse, int nv, const double* V, const int* vuse,
                  CondCovariates& cov) {
  if (nw < 0 || nv < 0 || K < 2) return -1;
  const std::size_t Kc = static_cast<std::size_t>(K - 1);
  cov = CondCovariates();
  cov.wcol.assign(Kc, {});
  cov.vcol.assign(Kc, {});
  if (nw > 0) {
    if (!W || !wuse) return -1;
    if (!finite_all(W, static_cast<std::size_t>(n_rows) * nw)) return -2;
    cov.W = W; cov.nw = static_cast<std::size_t>(nw);
    for (std::size_t k = 0; k < Kc; ++k)
      for (int j = 0; j < nw; ++j) {
        const int u = wuse[k * static_cast<std::size_t>(nw) + j];
        if (u != 0 && u != 1) return -1;
        if (u) cov.wcol[k].push_back(static_cast<std::size_t>(j));
      }
  }
  if (nv > 0) {
    if (!V || !vuse) return -1;
    if (!finite_all(V, static_cast<std::size_t>(n_rows) * nv)) return -2;
    cov.V = V; cov.nv = static_cast<std::size_t>(nv);
    for (std::size_t k = 0; k < Kc; ++k)
      for (int j = 0; j < nv; ++j) {
        const int u = vuse[k * static_cast<std::size_t>(nv) + j];
        if (u != 0 && u != 1) return -1;
        if (u) cov.vcol[k].push_back(static_cast<std::size_t>(j));
      }
  }
  return 0;
}
}  // namespace

extern "C" {

int sc_constrainedK_fit_cov(int n, int p, const double* time, const int* cause, const double* X, int K,
                            int alpha_free, double alpha_value, int baseline, const int* family,
                            int nw, const double* W, const int* wuse, int nv, const double* V, const int* vuse,
                            const double* init, int max_iter, double tol,
                            double* theta, double* var, sc_fit_info* info) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, 0, nullptr, s)) return rc;
  if (int rc = check_constrainedK(n, p, time, cause, X, K, family, s)) return rc;
  if (!theta) return -1;
  CondCovariates cov;
  if (int rc = make_cond_cov(n, K, nw, W, wuse, nv, V, vuse, cov)) return rc;
  const std::vector<int> fam(family, family + (K - 1));
  int ne = 0;
  for (int i = 0; i < n; ++i) ne += cause[i] != 0;
  if (cov.empty() || init) {                              // as given (start values are on the scale of W and V)
    ConstrainedKModel m(n, time, cause, X, s, K, fam, cov);
    const std::size_t P = m.nparam();
    std::vector<double> st;
    if (init) { if (!finite_all(init, P)) return -2; st.assign(init, init + P); }
    const FitResult f = constrainedK_fit(m, n, time, cause, fam, init ? &st : nullptr, max_iter, tol);
    export_fit(f, P, ne, theta, var, info);
    return 0;
  }
  // fitted with centred covariates (better conditioned), then expressed for W and V themselves
  std::vector<double> mw(cov.nw, 0.0), mv(cov.nv, 0.0), Wc, Vc;
  auto centre = [n](const double* A, std::size_t q, std::vector<double>& mean, std::vector<double>& out) {
    out.assign(A, A + static_cast<std::size_t>(n) * q);
    for (std::size_t j = 0; j < q; ++j) {
      double sm = 0.0;
      for (int i = 0; i < n; ++i) sm += A[static_cast<std::size_t>(i) * q + j];
      mean[j] = sm / n;
      for (int i = 0; i < n; ++i) out[static_cast<std::size_t>(i) * q + j] -= mean[j];
    }
  };
  if (cov.W) { centre(cov.W, cov.nw, mw, Wc); cov.W = Wc.data(); }
  if (cov.V) { centre(cov.V, cov.nv, mv, Vc); cov.V = Vc.data(); }
  ConstrainedKModel m(n, time, cause, X, s, K, fam, cov);
  const std::size_t P = m.nparam();
  FitResult f = constrainedK_fit(m, n, time, cause, fam, nullptr, max_iter, tol);
  if (f.beta.size() == P) {
    std::vector<double> J;
    cond_uncentre(m, mw, mv, f.beta, J);
    if (f.var.n == P) {
      Mat V2(P);
      for (std::size_t a = 0; a < P; ++a)
        for (std::size_t b = 0; b < P; ++b) {
          double v = 0.0;
          for (std::size_t c = 0; c < P; ++c) {
            if (J[a * P + c] == 0.0) continue;
            for (std::size_t d = 0; d < P; ++d) if (J[b * P + d] != 0.0) v += J[a * P + c] * f.var(c, d) * J[b * P + d];
          }
          V2(a, b) = v;
        }
      f.var = V2;
    }
  }
  export_fit(f, P, ne, theta, var, info);
  return 0;
} catch (...) { return -9; }

int sc_constrainedK_fit(int n, int p, const double* time, const int* cause, const double* X, int K,
                        int alpha_free, double alpha_value, int baseline, const int* family,
                        const double* init, int max_iter, double tol,
                        double* theta, double* var, sc_fit_info* info) {
  return sc_constrainedK_fit_cov(n, p, time, cause, X, K, alpha_free, alpha_value, baseline, family, 0, nullptr, nullptr,
                                 0, nullptr, nullptr, init, max_iter, tol, theta, var, info);
}

int sc_constrainedK_eval_cov(int n, int p, const double* time, const int* cause, const double* X, int K,
                             int alpha_free, double alpha_value, int baseline, const int* family,
                             int nw, const double* W, const int* wuse, int nv, const double* V, const int* vuse,
                             const double* theta, double* loglik, double* grad, double* info_matrix) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, 0, nullptr, s)) return rc;
  if (int rc = check_constrainedK(n, p, time, cause, X, K, family, s); rc != 0 && rc != -3) return rc;
  if (!theta) return -1;
  CondCovariates cov;
  if (int rc = make_cond_cov(n, K, nw, W, wuse, nv, V, vuse, cov)) return rc;
  ConstrainedKModel m(n, time, cause, X, s, K, std::vector<int>(family, family + (K - 1)), cov);
  Eval ev;
  m.eval(std::vector<double>(theta, theta + m.nparam()), ev, true);
  export_eval(ev, m.nparam(), loglik, grad, info_matrix);
  return 0;
} catch (...) { return -9; }

int sc_constrainedK_eval(int n, int p, const double* time, const int* cause, const double* X, int K,
                         int alpha_free, double alpha_value, int baseline, const int* family,
                         const double* theta, double* loglik, double* grad, double* info_matrix) {
  return sc_constrainedK_eval_cov(n, p, time, cause, X, K, alpha_free, alpha_value, baseline, family, 0, nullptr, nullptr,
                                  0, nullptr, nullptr, theta, loglik, grad, info_matrix);
}

int sc_constrainedK_predict_cov(int p, int K, int alpha_free, double alpha_value, int baseline, const int* family,
                                int nw, const int* wuse, int nv, const int* vuse,
                                const double* theta, const double* var, int k, int average, int n_new, const double* Xnew,
                                const double* Wnew, const double* Vnew,
                                int nt, const double* times, double* out_cif, double* out_se) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, 0, nullptr, s)) return rc;
  if (!theta || !family || !times || !out_cif || n_new <= 0 || nt <= 0 || K < 2 || k < 1 || k > K || (p > 0 && !Xnew)) return -1;
  for (int j = 0; j < K - 1; ++j) if (family[j] < 0 || family[j] >= COND_COUNT) return -1;
  if (out_se && !var) return -1;
  CondCovariates cov;
  if (int rc = make_cond_cov(n_new, K, nw, Wnew, wuse, nv, Vnew, vuse, cov)) return rc;
  std::vector<double> dummy_t(1, 0.0); std::vector<int> dummy_c(1, 0);
  ConstrainedKModel m(0, dummy_t.data(), dummy_c.data(), Xnew, s, K, std::vector<int>(family, family + (K - 1)), cov);
  const std::size_t P = m.nparam();
  const std::vector<double> th(theta, theta + P);
  const std::size_t nws = cov.W ? cov.nw : 0, nvs = cov.V ? cov.nv : 0;
  auto wr = [&](int i) { return cov.W ? cov.W + static_cast<std::size_t>(i) * nws : nullptr; };
  auto vr = [&](int i) { return cov.V ? cov.V + static_cast<std::size_t>(i) * nvs : nullptr; };
  std::vector<double> g;
  if (average) {                                         // mean over the rows; delta method for the mean
    std::vector<double> gs(P);
    for (int t = 0; t < nt; ++t) {
      double Fs = 0.0;
      std::fill(gs.begin(), gs.end(), 0.0);
      for (int i = 0; i < n_new; ++i) {
        Fs += m.cifk_grad(th, k, Xnew + static_cast<std::size_t>(i) * p, times[t], g, wr(i), vr(i));
        for (std::size_t a = 0; a < P; ++a) gs[a] += g[a];
      }
      out_cif[t] = Fs / n_new;
      if (out_se) {
        double v = 0.0;
        for (std::size_t a = 0; a < P; ++a) { if (gs[a] == 0.0) continue; for (std::size_t b = 0; b < P; ++b) v += gs[a] * var[a * P + b] * gs[b]; }
        v /= static_cast<double>(n_new) * n_new;
        out_se[t] = v >= 0.0 ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN();
      }
    }
    return 0;
  }
  for (int i = 0; i < n_new; ++i)
    for (int t = 0; t < nt; ++t) {
      const double F = m.cifk_grad(th, k, Xnew + static_cast<std::size_t>(i) * p, times[t], g, wr(i), vr(i));
      out_cif[static_cast<std::size_t>(i) * nt + t] = F;
      if (out_se) {
        double v = 0.0;
        for (std::size_t a = 0; a < P; ++a) { if (g[a] == 0.0) continue; for (std::size_t b = 0; b < P; ++b) v += g[a] * var[a * P + b] * g[b]; }
        out_se[static_cast<std::size_t>(i) * nt + t] = v >= 0.0 ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN();
      }
    }
  return 0;
} catch (...) { return -9; }

int sc_constrainedK_predict(int p, int K, int alpha_free, double alpha_value, int baseline, const int* family,
                            const double* theta, const double* var, int k, int average, int n_new, const double* Xnew,
                            int nt, const double* times, double* out_cif, double* out_se) {
  return sc_constrainedK_predict_cov(p, K, alpha_free, alpha_value, baseline, family, 0, nullptr, 0, nullptr, theta, var, k,
                                     average, n_new, Xnew, nullptr, nullptr, nt, times, out_cif, out_se);
}

int sc_subdist_fit(int n, int p, const double* time, const int* cause, const double* X, int interest,
                   int alpha_free, double alpha_value, int baseline, int n_knots, const double* knots, int grid,
                   const double* init, int max_iter, double tol,
                   double* theta, double* var_robust, double* var_model, sc_fit_info* info, int* n_support) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (int rc = check_subdist(n, p, time, cause, X, interest, s)) return rc;
  if (!theta || grid < 0) return -1;
  const std::size_t P = s.nparam();
  std::vector<int> c12(n);
  int ne = 0;
  for (int i = 0; i < n; ++i) { c12[i] = cause[i] == 0 ? 0 : (cause[i] == interest ? 1 : 2); ne += c12[i] == 1; }
  SubdistCIFModel m(n, time, c12.data(), X, s, 1, grid);
  std::vector<double> st;
  if (init) { if (!finite_all(init, P)) return -2; st.assign(init, init + P); }
  const SubdistFit f = subdist_fit(m, n, time, c12.data(), init ? &st : nullptr, max_iter, tol);
  export_fit(f.fit, P, ne, theta, var_robust, info);
  if (var_model) for (std::size_t j = 0; j < P * P; ++j) var_model[j] = f.var_model.n == P ? f.var_model.a[j] : std::numeric_limits<double>::quiet_NaN();
  if (n_support) *n_support = static_cast<int>(m.n_support());
  return 0;
} catch (...) { return -9; }

int sc_subdist_var_ipcw(int n, int p, const double* time, const int* cause, const double* X, int interest,
                        int alpha_free, double alpha_value, int baseline, int n_knots, const double* knots,
                        const double* theta, const double* var_model, double* var) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (int rc = check_subdist(n, p, time, cause, X, interest, s)) return rc;
  if (!theta || !var_model || !var) return -1;
  const std::size_t P = s.nparam();
  if (!finite_all(theta, P) || !finite_all(var_model, P * P)) return -2;
  std::vector<int> c12(n);
  for (int i = 0; i < n; ++i) c12[i] = cause[i] == 0 ? 0 : (cause[i] == interest ? 1 : 2);
  SubdistCIFModel m(n, time, c12.data(), X, s, 1, 0);
  Mat M;
  if (!m.meat_ipcw(std::vector<double>(theta, theta + P), M)) return -4;
  std::vector<double> AM(P * P, 0.0);
  for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) { double v = 0.0; for (std::size_t c = 0; c < P; ++c) v += var_model[a * P + c] * M(c, b); AM[a * P + b] = v; }
  for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) { double v = 0.0; for (std::size_t c = 0; c < P; ++c) v += AM[a * P + c] * var_model[c * P + b]; var[a * P + b] = v; }
  return 0;
} catch (...) { return -9; }

int sc_subdist_eval(int n, int p, const double* time, const int* cause, const double* X, int interest,
                    int alpha_free, double alpha_value, int baseline, int n_knots, const double* knots, int grid,
                    const double* theta, double* loglik, double* grad, double* info_matrix) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (int rc = check_subdist(n, p, time, cause, X, interest, s); rc != 0 && rc != -3) return rc;
  if (!theta || grid < 0) return -1;
  SubdistCIFModel m(n, time, cause, X, s, interest, grid);
  Eval ev;
  m.eval(std::vector<double>(theta, theta + s.nparam()), ev, true);
  export_eval(ev, s.nparam(), loglik, grad, info_matrix);
  return 0;
} catch (...) { return -9; }

int sc_constrained_fit(int n, int p, const double* time, const int* cause, const double* X,
                       int alpha_free, double alpha_value, int baseline, int family,
                       const double* init, int max_iter, double tol,
                       double* theta, double* var, sc_fit_info* info) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, 0, nullptr, s)) return rc;
  if (int rc = check_constrained(n, p, time, cause, X, family, s)) return rc;
  if (!theta) return -1;
  ConstrainedCIFModel m(n, time, cause, X, s, family);
  const std::size_t P = m.nparam();
  std::vector<double> st;
  if (init) { if (!finite_all(init, P)) return -2; st.assign(init, init + P); }
  const FitResult f = constrained_fit(m, n, time, cause, init ? &st : nullptr, max_iter, tol);
  int ne = 0;
  for (int i = 0; i < n; ++i) ne += cause[i] != 0;
  export_fit(f, P, ne, theta, var, info);
  return 0;
} catch (...) { return -9; }

int sc_constrained_eval(int n, int p, const double* time, const int* cause, const double* X,
                        int alpha_free, double alpha_value, int baseline, int family,
                        const double* theta, double* loglik, double* grad, double* info_matrix) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, 0, nullptr, s)) return rc;
  if (int rc = check_constrained(n, p, time, cause, X, family, s); rc != 0 && rc != -3) return rc;
  if (!theta) return -1;
  ConstrainedCIFModel m(n, time, cause, X, s, family);
  Eval ev;
  m.eval(std::vector<double>(theta, theta + m.nparam()), ev, true);
  export_eval(ev, m.nparam(), loglik, grad, info_matrix);
  return 0;
} catch (...) { return -9; }

int sc_constrained_predict(int p, int alpha_free, double alpha_value, int baseline, int family,
                           const double* theta, const double* var, int k, int n_new, const double* Xnew,
                           int nt, const double* times, double* out_cif, double* out_se) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, 0, nullptr, s)) return rc;
  if (!theta || !times || !out_cif || n_new <= 0 || nt <= 0 || k < 1 || k > 2 || (p > 0 && !Xnew)) return -1;
  if (family < 0 || family >= COND_COUNT || (out_se && !var)) return -1;
  std::vector<double> dummy_t(1, 0.0); std::vector<int> dummy_c(1, 0);
  ConstrainedCIFModel m(0, dummy_t.data(), dummy_c.data(), Xnew, s, family);
  const std::size_t P = m.nparam();
  const std::vector<double> th(theta, theta + P);
  std::vector<double> g;
  for (int i = 0; i < n_new; ++i)
    for (int t = 0; t < nt; ++t) {
      const double F = m.cif12_grad(th, k, Xnew + static_cast<std::size_t>(i) * p, times[t], g);
      out_cif[static_cast<std::size_t>(i) * nt + t] = F;
      if (out_se) {
        double v = 0.0;
        for (std::size_t a = 0; a < P; ++a) { if (g[a] == 0.0) continue; for (std::size_t b = 0; b < P; ++b) v += g[a] * var[a * P + b] * g[b]; }
        out_se[static_cast<std::size_t>(i) * nt + t] = v >= 0.0 ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN();
      }
    }
  return 0;
} catch (...) { return -9; }

}  // extern "C"

// ---------------------------------------------------------------------------------------------
// Phase 2c: curves for plotting -- averaged predictions and bootstrap intervals
// ---------------------------------------------------------------------------------------------
namespace {
// mean over the rows of Xnew of the CIF and of its gradient; delta-method standard error of the mean
template <class Model>
void mean_curve(const Model& m, int which, std::size_t P, const std::vector<double>& th, const double* var, int p, int n_new,
                const double* Xnew, int nt, const double* times, double* out_cif, double* out_se, bool two_cause) {
  std::vector<double> g, gs(P);
  for (int t = 0; t < nt; ++t) {
    double Fs = 0.0;
    std::fill(gs.begin(), gs.end(), 0.0);
    for (int i = 0; i < n_new; ++i) {
      const double* x = Xnew + static_cast<std::size_t>(i) * p;
      double F;
      if constexpr (std::is_same<Model, ConstrainedCIFModel>::value) F = m.cif12_grad(th, which, x, times[t], g);
      else F = m.cif_grad(th, which, x, times[t], g);
      (void)two_cause;
      Fs += F;
      for (std::size_t a = 0; a < P; ++a) gs[a] += g[a];
    }
    Fs /= n_new;
    for (double& v : gs) v /= n_new;
    out_cif[t] = Fs;
    if (out_se) {
      double v = 0.0;
      for (std::size_t a = 0; a < P; ++a) { if (gs[a] == 0.0) continue; for (std::size_t b = 0; b < P; ++b) v += gs[a] * var[a * P + b] * gs[b]; }
      out_se[t] = v >= 0.0 ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN();
    }
  }
}

double step_at(const std::vector<double>& x, const std::vector<double>& y, double at) {
  const std::size_t k = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), at) - x.begin());
  return k == 0 ? 0.0 : y[k - 1];
}

// CIF of `interest` at `times` for the rows of Xnew (or their mean when average) from one data set.
// model 0: Fine-Gray; 1: cause-specific Cox. adjust = 0: no covariates (all coefficients zero).
// Returns false when the model cannot be fitted to these data.
bool cif_once(int model, int n, int p, const double* time, const int* cause, const double* X, int K, int interest, int ties,
              int product_limit, int adjust, int n_new, const double* Xnew, int average, int nt, const double* times,
              std::vector<double>& out) {
  const int rows = average ? 1 : n_new;
  out.assign(static_cast<std::size_t>(rows) * nt, 0.0);
  sc_fit_info info;
  std::vector<double> beta(static_cast<std::size_t>(p) * (model == 1 ? K : 1), 0.0), v(static_cast<std::size_t>(p) * p);
  if (model == 0) {
    if (adjust) {
      std::vector<double> vr(static_cast<std::size_t>(p) * p);
      if (sc_finegray(n, p, time, cause, X, interest, 50, 1e-10, beta.data(), v.data(), vr.data(), &info) != 0) return false;
    }
    std::vector<double> bt(n), bh(n);
    int m = 0;
    if (sc_finegray_baseline(n, p, time, cause, X, interest, beta.data(), nullptr, bt.data(), bh.data(), &m) != 0) return false;
    bt.resize(m); bh.resize(m);
    for (int i = 0; i < n_new; ++i) {
      double lp = 0.0;
      for (int j = 0; j < p; ++j) lp += Xnew[static_cast<std::size_t>(i) * p + j] * beta[j];
      const double r = std::exp(lp);
      for (int t = 0; t < nt; ++t) {
        const double F = -std::expm1(-step_at(bt, bh, times[t]) * r);
        if (average) out[t] += F / n_new; else out[static_cast<std::size_t>(i) * nt + t] = F;
      }
    }
  } else {
    std::vector<int> st(n);
    if (adjust)
      for (int k = 0; k < K; ++k) {
        for (int i = 0; i < n; ++i) st[i] = cause[i] == k + 1;
        if (sc_coxph(n, p, time, st.data(), X, nullptr, ties, 50, 1e-10, beta.data() + static_cast<std::size_t>(k) * p, v.data(), &info) != 0)
          return false;
      }
    std::vector<double> grid(n), cif;
    cif.assign(static_cast<std::size_t>(n_new) * n * K, 0.0);
    int m = 0;
    if (sc_csc_cif(n, p, time, cause, X, K, beta.data(), ties, n_new, Xnew, product_limit, grid.data(), cif.data(), &m) != 0) return false;
    grid.resize(m);
    for (int t = 0; t < nt; ++t) {
      const std::size_t q = static_cast<std::size_t>(std::upper_bound(grid.begin(), grid.end(), times[t]) - grid.begin());
      for (int i = 0; i < n_new; ++i) {
        const double F = q == 0 ? 0.0 : cif[(static_cast<std::size_t>(i) * m + (q - 1)) * K + (interest - 1)];
        if (average) out[t] += F / n_new; else out[static_cast<std::size_t>(i) * nt + t] = F;
      }
    }
  }
  for (double vv : out) if (!std::isfinite(vv)) return false;
  return true;
}
}  // namespace

extern "C" {

int sc_directcif_predict_mean(int p, int K, const int* alpha_free, const double* alpha_value, const int* baseline,
                              const int* n_knots, const double* knots, const double* theta, const double* var, int k,
                              int n_new, const double* Xnew, int nt, const double* times, double* out_cif, double* out_se) try {
  DirectSpec s;
  if (int rc = make_spec(p, K, alpha_free, alpha_value, baseline, n_knots, knots, s)) return rc;
  if (!theta || !times || !out_cif || n_new <= 0 || nt <= 0 || k < 1 || k > K || (p > 0 && !Xnew) || (out_se && !var)) return -1;
  std::vector<double> dummy_t(1, 0.0); std::vector<int> dummy_c(1, 0);
  DirectCIFModel m(0, dummy_t.data(), dummy_c.data(), Xnew, s);
  const std::size_t P = s.nparam();
  mean_curve(m, k - 1, P, std::vector<double>(theta, theta + P), var, p, n_new, Xnew, nt, times, out_cif, out_se, false);
  return 0;
} catch (...) { return -9; }

int sc_constrained_predict_mean(int p, int alpha_free, double alpha_value, int baseline, int family, const double* theta,
                                const double* var, int k, int n_new, const double* Xnew, int nt, const double* times,
                                double* out_cif, double* out_se) try {
  DirectSpec s;
  if (int rc = one_spec(p, alpha_free, alpha_value, baseline, 0, nullptr, s)) return rc;
  if (!theta || !times || !out_cif || n_new <= 0 || nt <= 0 || k < 1 || k > 2 || (p > 0 && !Xnew) || (out_se && !var)) return -1;
  if (family < 0 || family >= COND_COUNT) return -1;
  std::vector<double> dummy_t(1, 0.0); std::vector<int> dummy_c(1, 0);
  ConstrainedCIFModel m(0, dummy_t.data(), dummy_c.data(), Xnew, s, family);
  const std::size_t P = m.nparam();
  mean_curve(m, k, P, std::vector<double>(theta, theta + P), var, p, n_new, Xnew, nt, times, out_cif, out_se, true);
  return 0;
} catch (...) { return -9; }

int sc_cif_bootstrap(int model, int n, int p, const double* time, const int* cause, const double* X, int K, int interest,
                     int ties, int product_limit, int adjust, int n_new, const double* Xnew, int average, int nt,
                     const double* times, int B, unsigned long long seed, double level,
                     double* out_est, double* out_se, double* out_lower, double* out_upper, int* n_ok) try {
  if (int rc = check_common(n, p, time, X)) return rc;
  if ((model != 0 && model != 1) || !cause || K < 1 || interest < 1 || interest > K || n_new <= 0 || nt <= 0 || !times ||
      !out_est || B < 0 || !(level > 0.0 && level < 1.0) || (p > 0 && !Xnew)) return -1;
  if (adjust && p == 0) return -1;
  for (int i = 0; i < n; ++i) if (cause[i] < 0 || cause[i] > K) return -2;
  // without covariates the models are fitted with one column of zeros and coefficient zero
  const int pe = adjust ? p : 1;
  std::vector<double> X0, Xn0;
  const double* Xu = X; const double* Xnu = Xnew;
  if (!adjust) { X0.assign(n, 0.0); Xn0.assign(n_new, 0.0); Xu = X0.data(); Xnu = Xn0.data(); }
  std::vector<double> est;
  if (!cif_once(model, n, pe, time, cause, Xu, K, interest, ties, product_limit, adjust, n_new, Xnu, average, nt, times, est)) return -4;
  const std::size_t cells = est.size();
  for (std::size_t j = 0; j < cells; ++j) out_est[j] = est[j];
  if (n_ok) *n_ok = 0;
  if (B == 0) return 0;
  std::vector<double> rep;
  std::vector<std::vector<double>> draws(cells);
  std::vector<double> tb(n), Xb(static_cast<std::size_t>(n) * pe);
  std::vector<int> cb(n);
  std::uint64_t s = seed;
  int ok = 0;
  for (int b = 0; b < B; ++b) {
    for (int i = 0; i < n; ++i) {
      const std::size_t j = draw_index(s, static_cast<std::size_t>(n));
      tb[i] = time[j]; cb[i] = cause[j];
      for (int q = 0; q < pe; ++q) Xb[static_cast<std::size_t>(i) * pe + q] = Xu[j * pe + q];
    }
    if (!cif_once(model, n, pe, tb.data(), cb.data(), Xb.data(), K, interest, ties, product_limit, adjust, n_new, Xnu, average, nt, times, rep))
      continue;
    ++ok;
    for (std::size_t j = 0; j < cells; ++j) draws[j].push_back(rep[j]);
  }
  if (n_ok) *n_ok = ok;
  const double a = (1.0 - level) / 2.0;
  auto quant = [](const std::vector<double>& v, double pr) {          // type 7 of Hyndman & Fan (1996)
    const double h = (v.size() - 1) * pr;
    const std::size_t lo = static_cast<std::size_t>(std::floor(h));
    const std::size_t hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (h - lo) * (v[hi] - v[lo]);
  };
  for (std::size_t j = 0; j < cells; ++j) {
    std::vector<double>& v = draws[j];
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (v.size() < 2) { if (out_se) out_se[j] = nan; if (out_lower) out_lower[j] = nan; if (out_upper) out_upper[j] = nan; continue; }
    double m1 = 0.0, m2 = 0.0;
    for (double x : v) m1 += x;
    m1 /= v.size();
    for (double x : v) m2 += (x - m1) * (x - m1);
    if (out_se) out_se[j] = std::sqrt(m2 / (v.size() - 1));
    std::sort(v.begin(), v.end());
    if (out_lower) out_lower[j] = quant(v, a);
    if (out_upper) out_upper[j] = quant(v, 1.0 - a);
  }
  return 0;
} catch (...) { return -9; }

}  // extern "C"
