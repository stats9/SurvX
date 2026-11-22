// survcore: shared gamma frailty models (clustered survival data).
//
// Subjects j in cluster c share an unobserved frailty z_c ~ Gamma(mean 1, variance theta) that
// multiplies their hazards: h_j(t | z_c) = z_c h_0(t) exp(x_j' beta). Integrating z_c out gives the
// marginal likelihood of cluster c (Clayton 1978; Klein 1992, Biometrics 48:795-806):
//   L_c = prod_j h_j(t_j)^{d_j} * prod_{m=0}^{D_c-1} (1 + m theta) * (1 + theta Lambda_c)^{-(1/theta + D_c)},
// with D_c = sum_j d_j and Lambda_c = sum_j H_0(t_j) exp(x_j' beta). The posterior of z_c is
// Gamma(shape 1/theta + D_c, rate 1/theta + Lambda_c) (empirical Bayes frailty).
//
// Writing G_c(Lambda, theta) = sum_m log(1 + m theta) - (1/theta + D_c) log(1 + theta Lambda), its
// derivatives are closed form (below); theta enters as phi = log theta. Two baselines:
//  * parametric (ParFrailtyModel): any family of parametric.hpp; full Newton with exact derivatives;
//  * semi-parametric (cox_frailty_fit): Breslow baseline with jumps h_k at the distinct event times,
//    profiled out by the EM fixed point h_k = d_k / sum_{t_j >= t_k} E[z_c] exp(x_j' beta)
//    (Nielsen, Gill, Andersen & Sorensen 1992, Scand J Statist 19:25-43; Klein 1992); the profile
//    log-likelihood in (beta, phi) has the exact gradient of the full likelihood at the profiled
//    baseline (envelope theorem) and is maximised by Newton with the Hessian from central
//    differences of that gradient; its inverse is the variance (profile likelihood, Murphy & van der
//    Vaart 2000, JASA 95:449-465), which accounts for the estimation of theta.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "linalg.hpp"
#include "newton.hpp"
#include "parametric.hpp"

namespace survcore {

// derivatives of G_c with respect to Lambda and phi = log theta
struct FrailtyG {
  double G = 0, GL = 0, GLL = 0, Gp = 0, Gpp = 0, GLp = 0;
};

// g(x) = log(1 + x)/x^2 - 1/(x (1 + x)) and its derivative, by series near 0
inline void frailty_gfun(double x, double& g, double& dg) {
  if (x < 1e-3) {
    g = 0.5 - 2.0 * x / 3.0 + 0.75 * x * x - 0.8 * x * x * x;
    dg = -2.0 / 3.0 + 1.5 * x - 2.4 * x * x;
  } else {
    const double l = std::log1p(x);
    g = l / (x * x) - 1.0 / (x * (1.0 + x));
    dg = 1.0 / ((1.0 + x) * x * x) - 2.0 * l / (x * x * x) + (1.0 + 2.0 * x) / (x * x * (1.0 + x) * (1.0 + x));
  }
}

inline FrailtyG frailty_G(double Lambda, int D, double theta) {
  FrailtyG r;
  const double x = theta * Lambda, opx = 1.0 + x;
  double s1 = 0, s2 = 0, s0 = 0;
  for (int m = 1; m < D; ++m) { const double a = 1.0 + m * theta; s0 += std::log(a); s1 += m / a; s2 += m * m / (a * a); }
  r.G = s0 - (1.0 / theta + D) * std::log1p(x);
  r.GL = -(1.0 + theta * D) / opx;
  r.GLL = theta * (1.0 + theta * D) / (opx * opx);
  double g, dg;
  frailty_gfun(x, g, dg);
  const double Gt = s1 + Lambda * Lambda * g - D * Lambda / opx;
  const double Gtt = -s2 + Lambda * Lambda * Lambda * dg + D * Lambda * Lambda / (opx * opx);
  const double GLt = -(D - Lambda) / (opx * opx);
  r.Gp = theta * Gt;
  r.Gpp = theta * theta * Gtt + theta * Gt;
  r.GLp = theta * GLt;
  return r;
}

// ------------------------------------------------------------------------------------------------
// parametric baseline; parameters (theta of the family, phi); cluster: 0-based ids, C clusters
// ------------------------------------------------------------------------------------------------
struct ParFrailtyModel {
  ParSpec sp;
  ParData dat;
  const int* cluster = nullptr;
  int C = 0;
  void eval(const std::vector<double>& th, Eval& ev, bool hess) const {
    const int P0 = sp.npar(), P = P0 + 1, p = sp.p, ns = sp.nshape(), q = 1 + ns;
    const double phi = th[P0], theta = std::exp(phi);
    const double* sh = th.data() + 1 + p;
    ev.ll = 0.0; ev.g.assign(P, 0.0);
    if (hess) ev.H = Mat(P);
    std::vector<double> Lam(C, 0.0), gL(static_cast<std::size_t>(C) * P0, 0.0), HL;
    std::vector<int> D(C, 0);
    if (hess) HL.assign(static_cast<std::size_t>(C) * P0 * P0, 0.0);
    std::vector<double> x(p + 1), g0(q), H0(q * q), g1(q), H1(q * q), gg(P0), hh(static_cast<std::size_t>(P0) * P0);
    double l0, l1;
    auto to_global = [&](const std::vector<double>& gl, const std::vector<double>& Hl, double sign) {
      std::fill(gg.begin(), gg.end(), 0.0);
      for (int a = 0; a <= p; ++a) gg[a] = sign * gl[0] * x[a];
      for (int c = 1; c < q; ++c) gg[p + c] = sign * gl[c];
      if (!hess) return;
      std::fill(hh.begin(), hh.end(), 0.0);
      for (int a = 0; a <= p; ++a) {
        for (int b = 0; b <= p; ++b) hh[a * P0 + b] = sign * Hl[0] * x[a] * x[b];
        for (int c = 1; c < q; ++c) { hh[(p + c) * P0 + a] = sign * Hl[c] * x[a]; hh[a * P0 + p + c] = hh[(p + c) * P0 + a]; }
      }
      for (int c = 1; c < q; ++c) for (int e = 1; e < q; ++e) hh[(p + c) * P0 + p + e] = sign * Hl[c * q + e];
    };
    for (int i = 0; i < dat.n; ++i) {
      const double* xr = dat.X + static_cast<std::size_t>(i) * p;
      x[0] = 1.0; for (int j = 0; j < p; ++j) x[j + 1] = xr[j];
      const double u = par_u(sp, xr, th);
      const int c = cluster[i];
      // H_i = -log S(t_i): from the contribution without an event
      par_local_all(sp, 0.0, dat.time[i], 0, u, sh, l0, g0.data(), H0.data());
      if (!std::isfinite(l0)) { ev.ll = -INFINITY; return; }
      Lam[c] += -l0;
      to_global(g0, H0, -1.0);
      for (int a = 0; a < P0; ++a) gL[static_cast<std::size_t>(c) * P0 + a] += gg[a];
      if (hess) for (int a = 0; a < P0 * P0; ++a) HL[static_cast<std::size_t>(c) * P0 * P0 + a] += hh[a];
      if (dat.status[i]) {
        ++D[c];
        // log h_i = l1 - l0
        par_local_all(sp, 0.0, dat.time[i], 1, u, sh, l1, g1.data(), H1.data());
        if (!std::isfinite(l1)) { ev.ll = -INFINITY; return; }
        ev.ll += l1 - l0;
        for (int k = 0; k < q; ++k) g1[k] -= g0[k];
        for (int k = 0; k < q * q; ++k) H1[k] -= H0[k];
        to_global(g1, H1, 1.0);
        for (int a = 0; a < P0; ++a) ev.g[a] += gg[a];
        if (hess) for (int a = 0; a < P0; ++a) for (int b = 0; b <= a; ++b) ev.H(a, b) -= hh[a * P0 + b];
      }
    }
    for (int c = 0; c < C; ++c) {
      const FrailtyG f = frailty_G(Lam[c], D[c], theta);
      if (!std::isfinite(f.G)) { ev.ll = -INFINITY; return; }
      ev.ll += f.G;
      const double* gc = &gL[static_cast<std::size_t>(c) * P0];
      for (int a = 0; a < P0; ++a) ev.g[a] += f.GL * gc[a];
      ev.g[P0] += f.Gp;
      if (hess) {
        const double* Hc = &HL[static_cast<std::size_t>(c) * P0 * P0];
        for (int a = 0; a < P0; ++a) {
          for (int b = 0; b <= a; ++b) ev.H(a, b) -= f.GLL * gc[a] * gc[b] + f.GL * Hc[a * P0 + b];
          ev.H(P0, a) -= f.GLp * gc[a];
        }
        ev.H(P0, P0) -= f.Gpp;
      }
    }
    if (!std::isfinite(ev.ll)) ev.ll = -INFINITY;
    if (hess) for (int a = 0; a < P; ++a) for (int b = a + 1; b < P; ++b) ev.H(a, b) = ev.H(b, a);
  }
  // posterior Gamma(shape, rate) of every cluster at th
  void posterior(const std::vector<double>& th, std::vector<double>& shape, std::vector<double>& rate) const {
    const int P0 = sp.npar(), p = sp.p, q = 1 + sp.nshape();
    const double theta = std::exp(th[P0]);
    const double* sh = th.data() + 1 + p;
    std::vector<double> Lam(C, 0.0), g0(q), H0(q * q);
    std::vector<int> D(C, 0);
    double l0;
    for (int i = 0; i < dat.n; ++i) {
      const double* xr = dat.X + static_cast<std::size_t>(i) * p;
      par_local_all(sp, 0.0, dat.time[i], 0, par_u(sp, xr, th), sh, l0, g0.data(), H0.data());
      Lam[cluster[i]] += -l0;
      D[cluster[i]] += dat.status[i] ? 1 : 0;
    }
    shape.resize(C); rate.resize(C);
    for (int c = 0; c < C; ++c) { shape[c] = 1.0 / theta + D[c]; rate[c] = 1.0 / theta + Lam[c]; }
  }
};

// ------------------------------------------------------------------------------------------------
// semi-parametric (Breslow) baseline
// ------------------------------------------------------------------------------------------------
struct CoxFrailtyData {
  int n = 0, p = 0, C = 0, K = 0;
  std::vector<double> time, X;             // sorted by time (ascending), X row-major
  std::vector<int> status, cluster;
  std::vector<double> tk;                  // distinct event times
  std::vector<int> dk;                     // events at each
  std::vector<int> idx;                    // number of distinct event times <= time_j
  std::vector<int> first;                  // first sorted row with time >= tk[k]
  std::vector<int> D;                      // events per cluster
};

inline CoxFrailtyData cox_frailty_data(int n, int p, const double* time, const int* status, const double* X, const int* cluster, int C) {
  CoxFrailtyData d;
  d.n = n; d.p = p; d.C = C;
  std::vector<int> o(n);
  for (int i = 0; i < n; ++i) o[i] = i;
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return time[a] < time[b]; });
  d.time.resize(n); d.status.resize(n); d.cluster.resize(n); d.X.resize(static_cast<std::size_t>(n) * p);
  for (int r = 0; r < n; ++r) {
    const int i = o[r];
    d.time[r] = time[i]; d.status[r] = status[i] ? 1 : 0; d.cluster[r] = cluster[i];
    for (int j = 0; j < p; ++j) d.X[static_cast<std::size_t>(r) * p + j] = X[static_cast<std::size_t>(i) * p + j];
  }
  for (int r = 0; r < n; ++r) {
    if (!d.status[r]) continue;
    if (d.tk.empty() || d.time[r] > d.tk.back()) { d.tk.push_back(d.time[r]); d.dk.push_back(1); }
    else d.dk.back() += 1;
  }
  d.K = static_cast<int>(d.tk.size());
  d.idx.resize(n);
  int k = 0;
  for (int r = 0; r < n; ++r) { while (k < d.K && d.tk[k] <= d.time[r]) ++k; d.idx[r] = k; }
  d.first.resize(d.K);
  int r = 0;
  for (int kk = 0; kk < d.K; ++kk) { while (r < n && d.time[r] < d.tk[kk]) ++r; d.first[kk] = r; }
  d.D.assign(C, 0);
  for (int rr = 0; rr < n; ++rr) d.D[d.cluster[rr]] += d.status[rr];
  return d;
}

struct CoxFrailtyState {
  std::vector<double> h;        // baseline hazard jumps
  std::vector<double> Lam, a;   // per cluster: Lambda_c and E[z_c]
  std::vector<double> e, H;     // per row: exp(x beta), H_0(t)
  int em_iter = 0;
};

// profile out the baseline at (beta, theta) by the EM fixed point; theta = 0: the Cox model
inline bool cox_frailty_profile(const CoxFrailtyData& d, const double* beta, double theta, CoxFrailtyState& s,
                                int max_em = 20000, double tol = 1e-13) {
  const int n = d.n, p = d.p, K = d.K, C = d.C;
  s.e.resize(n); s.H.resize(n); s.Lam.assign(C, 0.0); s.a.assign(C, 1.0);
  for (int r = 0; r < n; ++r) {
    double eta = 0;
    for (int j = 0; j < p; ++j) eta += d.X[static_cast<std::size_t>(r) * p + j] * beta[j];
    s.e[r] = std::exp(eta);
    if (!std::isfinite(s.e[r])) return false;
  }
  if (static_cast<int>(s.h.size()) != K) s.h.assign(K, 0.0);
  std::vector<double> cum(K + 1), w(n), hn(K);
  bool first_pass = s.h.empty() || s.h[0] <= 0.0;
  for (int it = 0; it < max_em; ++it) {
    s.em_iter = it + 1;
    if (!first_pass) {
      cum[0] = 0.0;
      for (int k = 0; k < K; ++k) cum[k + 1] = cum[k] + s.h[k];
      std::fill(s.Lam.begin(), s.Lam.end(), 0.0);
      for (int r = 0; r < n; ++r) { s.H[r] = cum[d.idx[r]]; s.Lam[d.cluster[r]] += s.e[r] * s.H[r]; }
      for (int c = 0; c < C; ++c) s.a[c] = theta > 0 ? (1.0 + theta * d.D[c]) / (1.0 + theta * s.Lam[c]) : 1.0;
    }
    // reverse cumulative risk sums
    double run = 0.0;
    std::vector<double> S(K);
    int r = n - 1;
    for (int k = K - 1; k >= 0; --k) {
      while (r >= d.first[k]) { run += s.a[d.cluster[r]] * s.e[r]; --r; }
      S[k] = run;
    }
    double change = 0.0;
    for (int k = 0; k < K; ++k) {
      hn[k] = d.dk[k] / S[k];
      if (!first_pass) change = std::fmax(change, std::fabs(hn[k] - s.h[k]) / hn[k]);
    }
    s.h = hn;
    if (theta <= 0.0 && !first_pass) break;          // Cox: the Breslow estimator in one step
    if (!first_pass && change < tol) break;
    first_pass = false;
  }
  // final H, Lambda, a at h
  cum[0] = 0.0;
  for (int k = 0; k < K; ++k) cum[k + 1] = cum[k] + s.h[k];
  std::fill(s.Lam.begin(), s.Lam.end(), 0.0);
  for (int r = 0; r < n; ++r) { s.H[r] = cum[d.idx[r]]; s.Lam[d.cluster[r]] += s.e[r] * s.H[r]; }
  for (int c = 0; c < C; ++c) s.a[c] = theta > 0 ? (1.0 + theta * d.D[c]) / (1.0 + theta * s.Lam[c]) : 1.0;
  return true;
}

// profile log-likelihood and its gradient in (beta, phi); theta = 0 (phi = -Inf): the Cox model
inline double cox_frailty_value(const CoxFrailtyData& d, const double* beta, double theta, CoxFrailtyState& s,
                                std::vector<double>* grad) {
  if (!cox_frailty_profile(d, beta, theta, s)) return -INFINITY;
  const int n = d.n, p = d.p, K = d.K, C = d.C;
  double ll = 0.0;
  for (int k = 0; k < K; ++k) ll += d.dk[k] * std::log(s.h[k]);
  for (int r = 0; r < n; ++r) if (d.status[r]) ll += std::log(s.e[r]);
  std::vector<double> GL(C, -1.0);
  double gphi = 0.0;
  for (int c = 0; c < C; ++c) {
    if (theta > 0) {
      const FrailtyG f = frailty_G(s.Lam[c], d.D[c], theta);
      ll += f.G; GL[c] = f.GL; gphi += f.Gp;
    } else {
      ll -= s.Lam[c];
    }
  }
  if (grad) {
    grad->assign(p + 1, 0.0);
    for (int r = 0; r < n; ++r) {
      const double wr = d.status[r] + GL[d.cluster[r]] * s.e[r] * s.H[r];
      for (int j = 0; j < p; ++j) (*grad)[j] += wr * d.X[static_cast<std::size_t>(r) * p + j];
    }
    (*grad)[p] = gphi;
  }
  return ll;
}

struct CoxFrailtyResult {
  std::vector<double> est;      // beta, phi
  Mat var;
  double ll = -INFINITY;
  int iterations = 0, flags = 0;
  bool converged = false;
  CoxFrailtyState state;
};

// Newton on the profile log-likelihood in (beta, phi), Hessian by central differences of the gradient
inline CoxFrailtyResult cox_frailty_fit(const CoxFrailtyData& d, std::vector<double> est, int max_iter, double tol) {
  CoxFrailtyResult r;
  const int p = d.p, P = p + 1;
  CoxFrailtyState s;
  std::vector<double> g, gp, gm;
  auto value = [&](const std::vector<double>& v, std::vector<double>* gr) {
    return cox_frailty_value(d, v.data(), std::exp(v[p]), s, gr);
  };
  auto hessian = [&](const std::vector<double>& v, Mat& Hm) -> bool {
    Hm = Mat(P);
    const std::vector<double> h_keep = s.h;
    for (int a = 0; a < P; ++a) {
      const double st = 1e-4 * std::fmax(1.0, std::fabs(v[a]));
      std::vector<double> vp = v, vm = v;
      vp[a] += st; vm[a] -= st;
      s.h = h_keep; if (!std::isfinite(value(vp, &gp))) return false;
      s.h = h_keep; if (!std::isfinite(value(vm, &gm))) return false;
      for (int b = 0; b < P; ++b) Hm(a, b) = -(gp[b] - gm[b]) / (2.0 * st);
    }
    s.h = h_keep;
    for (int a = 0; a < P; ++a) for (int b = 0; b < a; ++b) { const double m = 0.5 * (Hm(a, b) + Hm(b, a)); Hm(a, b) = m; Hm(b, a) = m; }
    return true;
  };
  double ll = value(est, &g);
  if (!std::isfinite(ll)) { r.flags |= FIT_NONFINITE; r.est = est; return r; }
  Mat Hm;
  for (int it = 1; it <= max_iter; ++it) {
    r.iterations = it;
    if (!hessian(est, Hm)) { r.flags |= FIT_NONFINITE; break; }
    ll = value(est, &g);
    std::vector<double> step;
    bool ok = false;
    double lambda = 0.0;
    for (int k = 0; k < 30 && !ok; ++k) {
      Mat L = Hm;
      if (lambda > 0) for (int j = 0; j < P; ++j) L(j, j) += lambda * std::fmax(std::fabs(Hm(j, j)), 1.0);
      if (cholesky(L)) { step = g; chol_solve(L, step); ok = true; }
      else lambda = lambda == 0.0 ? 1e-6 : lambda * 10.0;
    }
    if (!ok) { r.flags |= FIT_INFO_NOT_PD; break; }
    // keep phi steps moderate
    if (std::fabs(step[p]) > 2.0) { const double f = 2.0 / std::fabs(step[p]); for (double& v : step) v *= f; }
    double scale = 1.0, lt = -INFINITY;
    std::vector<double> trial(P);
    const std::vector<double> h_keep = s.h;
    bool acc = false;
    for (int hh = 0; hh < 40; ++hh) {
      for (int j = 0; j < P; ++j) trial[j] = est[j] + scale * step[j];
      s.h = h_keep;
      lt = value(trial, nullptr);
      if (std::isfinite(lt) && lt >= ll - 1e-12 * (std::fabs(ll) + 1.0)) { acc = true; break; }
      scale *= 0.5;
    }
    if (!acc) { s.h = h_keep; r.flags |= FIT_NONFINITE; break; }
    est = trial;
    const double dll = lt - ll;
    ll = lt;
    if (std::fabs(dll) <= tol * (std::fabs(ll) + tol) || est[p] < std::log(1e-8)) { r.converged = true; break; }
  }
  r.ll = value(est, &g);
  if (hessian(est, Hm)) {
    Mat L = Hm;
    if (cholesky(L)) r.var = chol_inverse(L); else r.flags |= FIT_INFO_NOT_PD;
  }
  value(est, &g);
  r.est = est;
  r.state = s;
  if (!r.converged && r.flags == FIT_OK) r.flags |= FIT_MAX_ITER;
  return r;
}

}  // namespace survcore
