// survcore: parametric cure models.
//   Mixture cure model (Boag 1949; Berkson & Gage 1952; Farewell 1982, Biometrics 38:1041-1046):
//        S(t | x, z) = pi(z) + {1 - pi(z)} S_u(t | x),
//   where pi(z) = expit(gamma'z) is the cured fraction and S_u the survival of the uncured, a
//   parametric model of parametric.hpp (any distribution, PH or AFT) with covariates x.
//   Non-mixture (promotion time, bounded cumulative hazard) cure model (Yakovlev & Tsodikov 1996;
//   Chen, Ibrahim & Sinha 1999, JASA 94:909-919):
//        S(t | x, z) = pi(z)^{F_u(t | x)} = exp{-theta(z) F_u(t | x)},  theta = -log pi,
//   with the same parametrisation of pi, so that pi(z) is again the cured fraction (as in
//   flexsurvcure, Amdahl 2022).
// Delayed entry: contributions are divided by S(entry | x, z) (Klein & Moeschberger 2003, Sec. 3.5).
// The latency enters through log f_u(t), log S_u(t) and log S_u(entry) and their exact local
// derivatives from par_local_all; the cure part is differentiated analytically.
#pragma once
#include <cmath>
#include <vector>
#include "parametric.hpp"

namespace survcore {

struct CureData {
  int n = 0;
  const double* entry = nullptr;
  const double* time = nullptr;
  const int* status = nullptr;
  const double* Xl = nullptr;   // n x pl, latency covariates
  const double* Xc = nullptr;   // n x pc, cure covariates
};

struct CureModel {
  ParSpec sp;        // latency; sp.p = pl
  int pc = 0;        // cure covariates (without the intercept)
  bool mixture = true;
  CureData dat;

  int npar() const { return sp.npar() + 1 + pc; }

  // value, gradient (q + 1) and Hessian ((q + 1)^2) of one contribution in the local variables
  // (u, shape_1..shape_ns, c), with c = gamma'z the logit of the cured fraction
  void local(double e, double t, int d, double u, const double* sh, double c,
             double& l, std::vector<double>& g, std::vector<double>& H) const {
    const int ns = sp.nshape(), q = 1 + ns, Q = q + 1;
    g.assign(Q, 0.0); H.assign(Q * Q, 0.0); l = 0.0;
    std::vector<double> gA(q), HA(q * q);
    const double pi = 1.0 / (1.0 + std::exp(-c));
    const double om = 1.0 - pi;                         // uncured probability
    const double pi1 = pi * om, pi2 = pi1 * (1.0 - 2.0 * pi);
    // adds s * (term in c and the latency component (L, gA, HA)) for the survival part at a time
    auto surv_part = [&](double L, double s) {
      const double S = std::exp(L);
      if (mixture) {
        const double D = pi + om * S;
        const double r = om * S / D;
        l += s * std::log(D);
        g[q] += s * pi1 * (1.0 - S) / D;
        H[q * Q + q] += s * (pi2 * (1.0 - S) * D - std::pow(pi1 * (1.0 - S), 2)) / (D * D);
        for (int a = 0; a < q; ++a) {
          g[a] += s * r * gA[a];
          const double cv = -pi1 * S * gA[a] / (D * D);
          H[a * Q + q] += s * cv; H[q * Q + a] += s * cv;
          for (int b = 0; b < q; ++b) H[a * Q + b] += s * (r * (gA[a] * gA[b] + HA[a * q + b]) - r * r * gA[a] * gA[b]);
        }
      } else {
        const double th = std::log1p(std::exp(-c));    // theta = -log pi
        l += s * (-th * (1.0 - S));
        g[q] += s * om * (1.0 - S);
        H[q * Q + q] += s * (-pi1 * (1.0 - S));
        for (int a = 0; a < q; ++a) {
          g[a] += s * th * S * gA[a];
          const double cv = -om * S * gA[a];
          H[a * Q + q] += s * cv; H[q * Q + a] += s * cv;
          for (int b = 0; b < q; ++b) H[a * Q + b] += s * th * S * (gA[a] * gA[b] + HA[a * q + b]);
        }
      }
    };
    double L;
    // survival part at t (censored), or the density at t (event)
    if (d) {
      par_local_all(sp, 0.0, t, 1, u, sh, L, gA.data(), HA.data());    // log f_u(t)
      if (!std::isfinite(L)) { l = -INFINITY; return; }
      l += L;
      for (int a = 0; a < q; ++a) { g[a] += gA[a]; for (int b = 0; b < q; ++b) H[a * Q + b] += HA[a * q + b]; }
      if (mixture) {
        l += std::log(om); g[q] += -pi; H[q * Q + q] += -pi1;
      } else {
        const double th = std::log1p(std::exp(-c)), th1 = -om, th2 = pi1;
        l += std::log(th); g[q] += th1 / th; H[q * Q + q] += (th2 * th - th1 * th1) / (th * th);
        par_local_all(sp, 0.0, t, 0, u, sh, L, gA.data(), HA.data());  // log S_u(t)
        surv_part(L, 1.0);                                             // -theta F_u(t)
      }
    } else {
      par_local_all(sp, 0.0, t, 0, u, sh, L, gA.data(), HA.data());
      surv_part(L, 1.0);
    }
    if (e > 0.0) {
      par_local_all(sp, 0.0, e, 0, u, sh, L, gA.data(), HA.data());    // log S_u(entry)
      surv_part(L, -1.0);
    }
  }

  void eval(const std::vector<double>& th, Eval& ev, bool hess) const {
    const int P = npar(), pl = sp.p, ns = sp.nshape(), q = 1 + ns, Q = q + 1;
    const int c0 = sp.npar();                            // index of the cure intercept
    const double* sh = th.data() + 1 + pl;
    ev.ll = 0.0; ev.g.assign(P, 0.0);
    if (hess) ev.H = Mat(P);
    std::vector<double> g, H;
    for (int i = 0; i < dat.n; ++i) {
      const double* xl = dat.Xl + static_cast<std::size_t>(i) * pl;
      const double* xc = dat.Xc + static_cast<std::size_t>(i) * pc;
      const double u = par_u(sp, xl, th);
      double c = th[c0];
      for (int j = 0; j < pc; ++j) c += xc[j] * th[c0 + 1 + j];
      double l;
      local(dat.entry ? dat.entry[i] : 0.0, dat.time[i], dat.status[i], u, sh, c, l, g, H);
      ev.ll += l;
      if (!std::isfinite(l)) continue;
      // theta index sets of each local variable: u -> (0, 1..pl) with x; shapes -> fixed; c -> (c0, ...) with z
      auto each = [&](int a, auto&& f) {
        if (a == 0) { f(0, 1.0); for (int j = 0; j < pl; ++j) f(1 + j, xl[j]); }
        else if (a < q) f(pl + a, 1.0);
        else { f(c0, 1.0); for (int j = 0; j < pc; ++j) f(c0 + 1 + j, xc[j]); }
      };
      for (int a = 0; a < Q; ++a) each(a, [&](int k, double m) { ev.g[k] += g[a] * m; });
      if (hess)
        for (int a = 0; a < Q; ++a)
          for (int b = 0; b < Q; ++b) {
            const double h = H[a * Q + b];
            if (h == 0.0) continue;
            each(a, [&](int k1, double m1) { each(b, [&](int k2, double m2) { ev.H(k1, k2) -= h * m1 * m2; }); });
          }
    }
    if (!std::isfinite(ev.ll)) ev.ll = -std::numeric_limits<double>::infinity();
  }
};

// survival of the cure model at t for latency row xl and cure row xc
inline double cure_surv(const CureModel& m, const std::vector<double>& th, const double* xl, const double* xc, double t,
                        double* hazard = nullptr, double* surv_u = nullptr, double* cured = nullptr) {
  const int c0 = m.sp.npar();
  double c = th[c0];
  for (int j = 0; j < m.pc; ++j) c += xc[j] * th[c0 + 1 + j];
  const double pi = 1.0 / (1.0 + std::exp(-c));
  if (cured) *cured = pi;
  const ParQ q = t > 0.0 ? par_hazard(m.sp, th, xl, t) : ParQ{0.0, -INFINITY};
  const double Su = std::exp(-q.H), hu = std::exp(q.logh);
  if (surv_u) *surv_u = Su;
  double S;
  if (m.mixture) {
    S = pi + (1.0 - pi) * Su;
    if (hazard) *hazard = (1.0 - pi) * hu * Su / S;
  } else {
    const double theta = std::log1p(std::exp(-c));
    S = std::exp(-theta * (1.0 - Su));
    if (hazard) *hazard = theta * hu * Su;
  }
  return S;
}

}  // namespace survcore
