// survcore C API: parametric survival models, cumulative incidence from parametric cause-specific
// hazards, and the kernel-smoothed hazard. See survcore/parametric.hpp and survcore/hazard.hpp.
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include "survcore/survcore.h"
#include "survcore/hazard.hpp"
#include "survcore/parametric.hpp"
#include "survcore/cure.hpp"
#include "survcore/interval.hpp"

using namespace survcore;

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();

bool spec_of(int dist, int aft, int p, int nk, const double* knots, ParSpec& sp) {
  sp.dist = dist; sp.aft = aft; sp.p = p;
  sp.knots.assign(knots && nk > 0 ? knots : nullptr, knots && nk > 0 ? knots + nk : nullptr);
  return sp.valid();
}

// one predicted quantity of a single model at a covariate row: 0 survival, 1 hazard, 2 cumulative
// hazard, 3 quantile (t = probability), 4 restricted mean survival time (t = tau)
double par_quantity(const ParSpec& sp, const std::vector<double>& th, const double* x, double t, int type) {
  switch (type) {
    case 0: return std::exp(-par_hazard(sp, th, x, t).H);
    case 1: return std::exp(par_hazard(sp, th, x, t).logh);
    case 2: return par_hazard(sp, th, x, t).H;
    case 3: return par_quantile(sp, th, x, t);
    default: return par_rmst(sp, th, x, t);
  }
}
}  // namespace

extern "C" {

int sc_par_fit(int n, int p, int dist, int aft, int nk, const double* knots, const double* entry, const double* time,
               const int* status, const double* X, const double* init, int max_iter, double tol, double* theta,
               double* var, sc_fit_info* info) try {
  ParSpec sp;
  if (n <= 0 || p < 0 || !time || !status || (p > 0 && !X) || !spec_of(dist, aft, p, nk, knots, sp)) return -1;
  int events = 0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(time[i]) || !(time[i] > 0.0)) return -2;
    if (entry && (!std::isfinite(entry[i]) || entry[i] < 0.0 || entry[i] >= time[i])) return -2;
    for (int j = 0; j < p; ++j) if (!std::isfinite(X[static_cast<std::size_t>(i) * p + j])) return -2;
    if (status[i]) ++events;
  }
  if (!events) return -3;
  ParData dat; dat.n = n; dat.entry = entry; dat.time = time; dat.status = status; dat.X = X;
  ParModel m{sp, dat};
  std::vector<double> start = init ? std::vector<double>(init, init + sp.npar()) : par_start(sp, dat);
  if (!init && sp.dist == PD_GENGAMMA) {
    // start from the Weibull AFT fit, which is the generalized gamma with Q = 1
    ParSpec w = sp; w.dist = PD_WEIBULL; w.knots.clear();
    ParModel mw{w, dat};
    FitResult rw = par_newton(mw, par_start(w, dat), max_iter, tol);
    if (std::isfinite(rw.ll)) { start.assign(rw.beta.begin(), rw.beta.end()); start.push_back(1.0); }
  }
  FitResult r = par_newton(m, start, max_iter, tol);
  const std::size_t P = sp.npar();
  for (std::size_t j = 0; j < P; ++j) theta[j] = r.beta[j];
  if (var) for (std::size_t j = 0; j < P * P; ++j) var[j] = r.var.n == P ? r.var.a[j] : kNaN;
  if (info) {
    info->loglik_init = r.ll_init; info->loglik = r.ll; info->grad_maxabs = r.grad_maxabs;
    info->iterations = r.iterations; info->converged = r.converged; info->flags = r.flags; info->n_events = events;
  }
  return 0;
} catch (...) { return -9; }

int sc_par_fit_ic(int n, int p, int dist, int aft, int nk, const double* knots, const double* time, const double* time2,
                  const int* type, const double* X, const double* init, int max_iter, double tol, double* theta,
                  double* var, sc_fit_info* info) try {
  ParSpec sp;
  if (n <= 0 || p < 0 || !time || !time2 || !type || !init || (p > 0 && !X) || !spec_of(dist, aft, p, nk, knots, sp)) return -1;
  int events = 0;
  for (int i = 0; i < n; ++i) {
    const int ty = type[i];
    if (ty < 0 || ty > 3) return -1;
    if (ty == 3) { if (!std::isfinite(time2[i]) || !(time2[i] > 0.0)) return -2; }
    else if (!std::isfinite(time[i]) || !(time[i] > 0.0)) return -2;
    if (ty == 2 && (!std::isfinite(time2[i]) || !(time2[i] > time[i]))) return -2;
    for (int j = 0; j < p; ++j) if (!std::isfinite(X[static_cast<std::size_t>(i) * p + j])) return -2;
    if (ty) ++events;
  }
  if (!events) return -3;
  ParData dat; dat.n = n; dat.time = time; dat.status = type; dat.X = X; dat.time2 = time2;
  ParModel m{sp, dat};
  FitResult r = par_newton(m, std::vector<double>(init, init + sp.npar()), max_iter, tol);
  const std::size_t P = sp.npar();
  for (std::size_t j = 0; j < P; ++j) theta[j] = r.beta[j];
  if (var) for (std::size_t j = 0; j < P * P; ++j) var[j] = r.var.n == P ? r.var.a[j] : kNaN;
  if (info) {
    info->loglik_init = r.ll_init; info->loglik = r.ll; info->grad_maxabs = r.grad_maxabs;
    info->iterations = r.iterations; info->converged = r.converged; info->flags = r.flags; info->n_events = events;
  }
  return 0;
} catch (...) { return -9; }

// Turnbull NPMLE: n_int receives the number of innermost intervals m (lo, hi, mass need room for n
// entries, m <= n); a, b the 0-based block of innermost intervals of every subject; w case weights
// (null: all one); info_d = (loglik, kkt), info_i = (iterations, converged)
int sc_turnbull(int n, const double* L, const double* R, const double* w, int max_iter, double tol, int* n_int,
                double* lo, double* hi, double* mass, int* a, int* b, double* info_d, int* info_i) try {
  if (n <= 0 || !L || !R || !n_int || !lo || !hi || !mass) return -1;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(L[i]) || L[i] < 0.0 || std::isnan(R[i]) || R[i] < L[i] || (R[i] == L[i] && !std::isfinite(R[i]))) return -2;
    if (w && (!std::isfinite(w[i]) || w[i] < 0.0)) return -2;
  }
  const TurnbullSupport s = turnbull_support(n, L, R);
  const TurnbullResult r = turnbull_fit(s, w, max_iter, tol);
  const int J = static_cast<int>(s.lo.size());
  if (J > n) return -9;
  *n_int = J;
  for (int j = 0; j < J; ++j) { lo[j] = s.lo[j]; hi[j] = s.hi[j]; mass[j] = r.p.empty() ? kNaN : r.p[j]; }
  if (a && b) for (int i = 0; i < n; ++i) { a[i] = s.a[i]; b[i] = s.b[i]; }
  if (info_d) { info_d[0] = r.loglik; info_d[1] = r.kkt; }
  if (info_i) { info_i[0] = r.iterations; info_i[1] = r.converged ? 1 : 0; }
  return 0;
} catch (...) { return -9; }

int sc_par_predict(int dist, int aft, int nk, const double* knots, int p, const double* theta, const double* var,
                   int n_new, const double* Xnew, int nt, const double* times, int type, int average,
                   double* out, double* out_se) try {
  ParSpec sp;
  if (!spec_of(dist, aft, p, nk, knots, sp) || n_new <= 0 || nt < 0 || type < 0 || type > 4) return -1;
  const std::size_t P = sp.npar();
  std::vector<double> th(theta, theta + P), V;
  if (var) V.assign(var, var + P * P);
  for (int j = 0; j < nt; ++j) {
    auto value = [&](const std::vector<double>& t, int i) { return par_quantity(sp, t, Xnew + static_cast<std::size_t>(i) * p, times[j], type); };
    if (average) {
      auto mean = [&](const std::vector<double>& t) { double s = 0; for (int i = 0; i < n_new; ++i) s += value(t, i); return s / n_new; };
      out[j] = mean(th);
      if (out_se) out_se[j] = var && std::isfinite(out[j]) ? delta_se(mean, th, V) : kNaN;
    } else {
      for (int i = 0; i < n_new; ++i) {
        const std::size_t c = static_cast<std::size_t>(i) * nt + j;
        out[c] = value(th, i);
        if (out_se) out_se[c] = var && std::isfinite(out[c]) ? delta_se([&](const std::vector<double>& t) { return value(t, i); }, th, V) : kNaN;
      }
    }
  }
  return 0;
} catch (...) { return -9; }

// K parametric cause-specific hazard models with the same covariates. dist, aft, nk: K codes and numbers of
// knots (cut points); knots: those of cause 1, then cause 2, ...; theta: the
// parameters of cause 1, then cause 2, ...; var: the K covariance matrices one after another (the
// models are fitted separately, so the joint covariance is block diagonal). k = 1..K: CIF of cause k;
// k = 0: overall survival. times must be increasing. out (and out_se): n_new x nt, or nt with average.
int sc_parcsc_predict(int K, const int* dist, const int* aft, const int* nk, const double* knots, int p,
                      const double* theta, const double* var, int k, int n_new, const double* Xnew, int nt,
                      const double* times, double tmax, int average, double* out, double* out_se) try {
  if (K < 1 || k < 0 || k > K || n_new <= 0 || nt <= 0) return -1;
  for (int j = 1; j < nt; ++j) if (!(times[j] >= times[j - 1])) return -1;
  std::vector<ParSpec> sp(K);
  std::vector<std::size_t> off(K + 1, 0);
  std::size_t koff = 0;
  for (int c = 0; c < K; ++c) {
    if (!spec_of(dist[c], aft[c], p, nk[c], knots ? knots + koff : nullptr, sp[c])) return -1;
    koff += nk[c] > 0 ? nk[c] : 0;
    off[c + 1] = off[c] + sp[c].npar();
  }
  const std::size_t Pt = off[K];
  std::vector<double> all(theta, theta + Pt);
  std::vector<double> tm(times, times + nt);
  for (double& t : tm) if (t < 0) t = 0;

  // the curve of every time for one row (or the average over rows) at parameters `a`
  auto curve = [&](const std::vector<double>& a, std::vector<double>& res) {
    std::vector<std::vector<double>> th(K);
    for (int c = 0; c < K; ++c) th[c].assign(a.begin() + off[c], a.begin() + off[c + 1]);
    res.assign(static_cast<std::size_t>(average ? 1 : n_new) * nt, 0.0);
    std::vector<double> cif(static_cast<std::size_t>(nt) * K), surv(nt);
    std::vector<const double*> xr(K);
    for (int i = 0; i < n_new; ++i) {
      for (int c = 0; c < K; ++c) xr[c] = Xnew + static_cast<std::size_t>(i) * p;
      par_cif(K, sp.data(), th.data(), xr.data(), nt, tm.data(), tmax, cif.data(), surv.data());
      for (int j = 0; j < nt; ++j) {
        const double v = k == 0 ? surv[j] : cif[static_cast<std::size_t>(j) * K + (k - 1)];
        if (average) res[j] += v / n_new; else res[static_cast<std::size_t>(i) * nt + j] = v;
      }
    }
  };
  std::vector<double> est;
  curve(all, est);
  std::copy(est.begin(), est.end(), out);
  if (out_se) {
    if (!var) { std::fill(out_se, out_se + est.size(), kNaN); return 0; }
    // central-difference gradient of every output in every parameter; then g' V g per block
    std::vector<std::vector<double>> G(Pt);
    std::vector<double> up, dn;
    for (std::size_t j = 0; j < Pt; ++j) {
      const double h = 1e-5 * std::fmax(1.0, std::fabs(all[j])), keep = all[j];
      all[j] = keep + h; curve(all, up);
      all[j] = keep - h; curve(all, dn);
      all[j] = keep;
      G[j].resize(est.size());
      for (std::size_t c = 0; c < est.size(); ++c) G[j][c] = (up[c] - dn[c]) / (2.0 * h);
    }
    std::size_t voff = 0;
    std::vector<double> v2(est.size(), 0.0);
    for (int c = 0; c < K; ++c) {
      const std::size_t Pc = sp[c].npar();
      for (std::size_t a = 0; a < Pc; ++a)
        for (std::size_t b = 0; b < Pc; ++b) {
          const double vab = var[voff + a * Pc + b];
          const std::vector<double>& ga = G[off[c] + a]; const std::vector<double>& gb = G[off[c] + b];
          for (std::size_t q = 0; q < est.size(); ++q) v2[q] += ga[q] * vab * gb[q];
        }
      voff += Pc * Pc;
    }
    for (std::size_t q = 0; q < est.size(); ++q) out_se[q] = v2[q] >= 0 && std::isfinite(v2[q]) ? std::sqrt(v2[q]) : kNaN;
  }
  return 0;
} catch (...) { return -9; }

// Kernel-smoothed hazard from the distinct event times et (increasing), events d and numbers at risk
// Y; bw <= 0 chooses the bandwidth by least-squares cross-validation (*bw_used receives it).
int sc_kernel_hazard(int m, const double* et, const double* d, const double* Y, double bw, int boundary,
                     int nx, const double* x, double* out_h, double* out_se, double* bw_used) try {
  if (m <= 0 || nx < 0 || !et || !d || !Y) return -1;
  for (int j = 0; j < m; ++j) if (!(Y[j] > 0) || !(d[j] >= 0) || (j && !(et[j] > et[j - 1]))) return -2;
  const double b = bw > 0 ? bw : hazard_bandwidth(m, et, d, Y, boundary != 0);
  if (bw_used) *bw_used = b;
  kernel_hazard(m, et, d, Y, b, boundary != 0, nx, x, out_h, out_se);
  return 0;
} catch (...) { return -9; }

// Parametric cure models (survcore/cure.hpp). theta = (latency: intercept, pl coefficients, shapes;
// cure: intercept, pc coefficients on the logit of the cured fraction).
int sc_cure_fit(int n, int pl, int pc, int dist, int aft, int nk, const double* knots, int mixture,
                const double* entry, const double* time, const int* status, const double* Xl, const double* Xc,
                const double* init, int max_iter, double tol, double* theta, double* var, sc_fit_info* info) try {
  ParSpec sp;
  if (n <= 0 || pl < 0 || pc < 0 || !time || !status || (pl > 0 && !Xl) || (pc > 0 && !Xc) || !spec_of(dist, aft, pl, nk, knots, sp)) return -1;
  int events = 0;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(time[i]) || !(time[i] > 0.0)) return -2;
    if (entry && (!std::isfinite(entry[i]) || entry[i] < 0.0 || entry[i] >= time[i])) return -2;
    for (int j = 0; j < pl; ++j) if (!std::isfinite(Xl[static_cast<std::size_t>(i) * pl + j])) return -2;
    for (int j = 0; j < pc; ++j) if (!std::isfinite(Xc[static_cast<std::size_t>(i) * pc + j])) return -2;
    if (status[i]) ++events;
  }
  if (!events) return -3;
  CureModel m;
  m.sp = sp; m.pc = pc; m.mixture = mixture != 0;
  m.dat.n = n; m.dat.entry = entry; m.dat.time = time; m.dat.status = status; m.dat.Xl = Xl; m.dat.Xc = Xc;
  const int P = m.npar();
  std::vector<double> start(P, 0.0);
  if (init) start.assign(init, init + P);
  else {
    // latency: the parametric model fitted to the event times alone (the uncured fail);
    // cure intercept: logit of the Kaplan-Meier estimate at the last time, kept within [0.02, 0.9]
    std::vector<int> idx;
    for (int i = 0; i < n; ++i) if (status[i]) idx.push_back(i);
    std::vector<double> te(idx.size()), xe(idx.size() * static_cast<std::size_t>(pl));
    std::vector<int> de(idx.size(), 1);
    for (std::size_t k = 0; k < idx.size(); ++k) {
      te[k] = time[idx[k]];
      for (int j = 0; j < pl; ++j) xe[k * pl + j] = Xl[static_cast<std::size_t>(idx[k]) * pl + j];
    }
    ParData de_dat; de_dat.n = static_cast<int>(idx.size()); de_dat.time = te.data(); de_dat.status = de.data(); de_dat.X = xe.data();
    ParSpec ws = sp;
    if (ws.dist == PD_GENGAMMA) ws.dist = PD_WEIBULL;
    ParModel lm{ws, de_dat};
    FitResult rl = par_newton(lm, par_start(ws, de_dat), max_iter, tol);
    std::vector<double> lat = std::isfinite(rl.ll) ? rl.beta : par_start(ws, de_dat);
    for (std::size_t j = 0; j < lat.size(); ++j) start[j] = lat[j];
    if (sp.dist == PD_GENGAMMA) start[1 + pl + 1] = 1.0;
    std::vector<int> o(n);
    std::iota(o.begin(), o.end(), 0);
    std::sort(o.begin(), o.end(), [&](int a, int b) { return time[a] < time[b]; });
    double S = 1.0; int risk = n;
    for (int k = 0; k < n;) {
      int j = k, d = 0;
      while (j < n && time[o[j]] == time[o[k]]) { d += status[o[j]] ? 1 : 0; ++j; }
      if (risk > 0) S *= 1.0 - static_cast<double>(d) / risk;
      risk -= j - k; k = j;
    }
    const double p0 = std::min(0.9, std::max(0.02, S));
    start[sp.npar()] = std::log(p0 / (1.0 - p0));
  }
  FitResult r = par_newton(m, start, max_iter, tol);
  for (int j = 0; j < P; ++j) theta[j] = r.beta[j];
  const std::size_t PP = static_cast<std::size_t>(P) * P;
  if (var) for (std::size_t j = 0; j < PP; ++j) var[j] = r.var.n == static_cast<std::size_t>(P) ? r.var.a[j] : kNaN;
  if (info) {
    info->loglik_init = r.ll_init; info->loglik = r.ll; info->grad_maxabs = r.grad_maxabs;
    info->iterations = r.iterations; info->converged = r.converged; info->flags = r.flags; info->n_events = events;
  }
  return 0;
} catch (...) { return -9; }

// type: 0 survival, 1 hazard, 2 survival of the uncured, 3 cured fraction, 4 cumulative hazard (-log S).
// out, out_se: n_new x nt (row-major), or nt with average.
int sc_cure_predict(int dist, int aft, int nk, const double* knots, int pl, int pc, int mixture, const double* theta,
                    const double* var, int n_new, const double* Xl, const double* Xc, int nt, const double* times,
                    int type, int average, double* out, double* out_se) try {
  CureModel m;
  if (!spec_of(dist, aft, pl, nk, knots, m.sp) || n_new <= 0 || nt <= 0 || type < 0 || type > 4) return -1;
  m.pc = pc; m.mixture = mixture != 0;
  const std::size_t P = m.npar();
  std::vector<double> th(theta, theta + P), V;
  if (var) V.assign(var, var + P * P);
  auto value = [&](const std::vector<double>& a, int i, double t) {
    const double* xl = Xl + static_cast<std::size_t>(i) * pl;
    const double* xc = Xc + static_cast<std::size_t>(i) * pc;
    double h, su, pi;
    const double S = cure_surv(m, a, xl, xc, t, &h, &su, &pi);
    switch (type) {
      case 0: return S;
      case 1: return h;
      case 2: return su;
      case 3: return pi;
      default: return -std::log(S);
    }
  };
  for (int j = 0; j < nt; ++j) {
    if (average) {
      auto mean = [&](const std::vector<double>& a) { double s = 0; for (int i = 0; i < n_new; ++i) s += value(a, i, times[j]); return s / n_new; };
      out[j] = mean(th);
      if (out_se) out_se[j] = var && std::isfinite(out[j]) ? delta_se(mean, th, V) : kNaN;
    } else {
      for (int i = 0; i < n_new; ++i) {
        const std::size_t c = static_cast<std::size_t>(i) * nt + j;
        out[c] = value(th, i, times[j]);
        if (out_se) out_se[c] = var && std::isfinite(out[c]) ? delta_se([&](const std::vector<double>& a) { return value(a, i, times[j]); }, th, V) : kNaN;
      }
    }
  }
  return 0;
} catch (...) { return -9; }

}  // extern "C"
