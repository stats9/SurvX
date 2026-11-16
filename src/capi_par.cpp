// survcore C API: parametric survival models, cumulative incidence from parametric cause-specific
// hazards, and the kernel-smoothed hazard. See survcore/parametric.hpp and survcore/hazard.hpp.
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include "survcore/survcore.h"
#include "survcore/hazard.hpp"
#include "survcore/parametric.hpp"

using namespace survcore;

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();

bool spec_of(int dist, int aft, int p, ParSpec& sp) {
  sp.dist = dist; sp.aft = aft; sp.p = p;
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

int sc_par_fit(int n, int p, int dist, int aft, const double* entry, const double* time, const int* status,
               const double* X, const double* init, int max_iter, double tol, double* theta, double* var,
               sc_fit_info* info) try {
  ParSpec sp;
  if (n <= 0 || p < 0 || !time || !status || (p > 0 && !X) || !spec_of(dist, aft, p, sp)) return -1;
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

int sc_par_predict(int dist, int aft, int p, const double* theta, const double* var, int n_new, const double* Xnew,
                   int nt, const double* times, int type, int average, double* out, double* out_se) try {
  ParSpec sp;
  if (!spec_of(dist, aft, p, sp) || n_new <= 0 || nt < 0 || type < 0 || type > 4) return -1;
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

// K parametric cause-specific hazard models with the same covariates. dist, aft: K codes; theta: the
// parameters of cause 1, then cause 2, ...; var: the K covariance matrices one after another (the
// models are fitted separately, so the joint covariance is block diagonal). k = 1..K: CIF of cause k;
// k = 0: overall survival. times must be increasing. out (and out_se): n_new x nt, or nt with average.
int sc_parcsc_predict(int K, const int* dist, const int* aft, int p, const double* theta, const double* var, int k,
                      int n_new, const double* Xnew, int nt, const double* times, double tmax, int average,
                      double* out, double* out_se) try {
  if (K < 1 || k < 0 || k > K || n_new <= 0 || nt <= 0) return -1;
  for (int j = 1; j < nt; ++j) if (!(times[j] >= times[j - 1])) return -1;
  std::vector<ParSpec> sp(K);
  std::vector<std::size_t> off(K + 1, 0);
  for (int c = 0; c < K; ++c) { if (!spec_of(dist[c], aft[c], p, sp[c])) return -1; off[c + 1] = off[c] + sp[c].npar(); }
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

}  // extern "C"
