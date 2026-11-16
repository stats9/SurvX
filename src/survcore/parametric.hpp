// survcore: parametric regression models for a single event time under right censoring and
// left truncation (delayed entry), and cumulative incidence functions from parametric
// cause-specific hazards.
//
// Families (Kalbfleisch & Prentice 2002, Ch. 2-3; Lawless 2003, Ch. 4-6):
//   proportional hazards (PH): H(t | z) = exp(eta) H0(t), eta = b0 + z'b
//     exponential  H0(t) = t
//     Weibull      H0(t) = t^kappa,              kappa = exp(s)        (flexsurv "weibullPH")
//     Gompertz     H0(t) = (exp(a t) - 1) / a,   a = s, any sign       (flexsurv "gompertz")
//   accelerated failure time (AFT): log T = mu + sigma W, mu = b0 + z'b, sigma = exp(s)
//     exponential  W standard minimum extreme value, sigma = 1         (survival::survreg)
//     Weibull      W standard minimum extreme value                     (survival::survreg)
//     lognormal    W standard normal                                    (survival::survreg)
//     loglogistic  W standard logistic                                  (survival::survreg)
// Parameter vector theta = (b0, b_1..b_p, s); s is absent for the exponential.
//
// Log-likelihood contribution of subject i with entry e_i >= 0, time t_i and status d_i:
//   d_i log h(t_i) + log S(t_i) - log S(e_i)                 (Klein & Moeschberger 2003, Sec. 3.5)
// Every contribution depends on theta only through two local quantities, u = eta or mu and s,
// so the gradient and Hessian are assembled exactly from the first and second partial derivatives
// in (u, s) by the chain rule; no numerical differentiation is used in the fit.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include "dual2.hpp"
#include "linalg.hpp"
#include "newton.hpp"

namespace survcore {

enum ParDist : int { PD_EXPONENTIAL = 0, PD_WEIBULL = 1, PD_GOMPERTZ = 2, PD_LOGNORMAL = 3, PD_LOGLOGISTIC = 4 };

struct ParSpec {
  int dist = PD_WEIBULL;
  int aft = 0;            // 0: PH, 1: AFT
  int p = 0;              // number of covariates (without the intercept)
  bool has_shape() const { return dist != PD_EXPONENTIAL; }
  int npar() const { return 1 + p + (has_shape() ? 1 : 0); }
  bool valid() const {
    if (dist < 0 || dist > 4 || p < 0) return false;
    if (aft == 0) return dist == PD_EXPONENTIAL || dist == PD_WEIBULL || dist == PD_GOMPERTZ;
    return dist == PD_EXPONENTIAL || dist == PD_WEIBULL || dist == PD_LOGNORMAL || dist == PD_LOGLOGISTIC;
  }
};

// value and partial derivatives of a contribution in the local quantities (u, s)
struct Loc { double l = 0, lu = 0, ls = 0, luu = 0, lus = 0, lss = 0; };

namespace par_detail {

constexpr double kLogSqrt2Pi = 0.91893853320467274178;

// log of the standard normal upper tail Phi(-w) and the inverse Mills ratio m(w) = phi(w)/Phi(-w)
inline void normal_upper(double w, double& logS, double& m) {
  if (w < 25.0) {
    const double S = 0.5 * std::erfc(w / std::sqrt(2.0));
    logS = std::log(S);
    m = std::exp(-0.5 * w * w - kLogSqrt2Pi) / S;
  } else {                                     // asymptotic expansion of the Mills ratio (Abramowitz & Stegun 26.2.12)
    const double w2 = w * w;
    const double r = 1.0 - 1.0 / w2 + 3.0 / (w2 * w2) - 15.0 / (w2 * w2 * w2);   // Phi(-w) w / phi(w)
    m = w / r;
    logS = -0.5 * w2 - kLogSqrt2Pi - std::log(w) + std::log(r);
  }
}

// log density g and log survivor G of the standardised error W of an AFT family, with derivatives
inline void aft_error(int dist, double w, double& g, double& g1, double& g2, double& G, double& G1, double& G2) {
  if (dist == PD_EXPONENTIAL || dist == PD_WEIBULL) {        // minimum extreme value
    const double e = std::exp(w);
    g = w - e; g1 = 1.0 - e; g2 = -e;
    G = -e; G1 = -e; G2 = -e;
  } else if (dist == PD_LOGNORMAL) {
    g = -0.5 * w * w - kLogSqrt2Pi; g1 = -w; g2 = -1.0;
    double logS, m; normal_upper(w, logS, m);
    G = logS; G1 = -m; G2 = -m * (m - w);
  } else {                                                    // logistic
    const double pr = w >= 0 ? 1.0 / (1.0 + std::exp(-w)) : std::exp(w) / (1.0 + std::exp(w));
    const double sp = w > 0 ? w + std::log1p(std::exp(-w)) : std::log1p(std::exp(w));   // log(1 + e^w)
    g = w - 2.0 * sp; g1 = 1.0 - 2.0 * pr; g2 = -2.0 * pr * (1.0 - pr);
    G = -sp; G1 = -pr; G2 = -pr * (1.0 - pr);
  }
}

// a term A(w), w = (y - u)/sigma, and its derivatives in (u, s)
inline void aft_term(double A, double A1, double A2, double w, double sigma, Loc& o, double sign) {
  o.l += sign * A;
  o.lu += sign * (-A1 / sigma);
  o.ls += sign * (-A1 * w);
  o.luu += sign * (A2 / (sigma * sigma));
  o.lus += sign * (A2 * w / sigma + A1 / sigma);
  o.lss += sign * (A2 * w * w + A1 * w);
}

// baseline cumulative hazard H0(t) of a PH family with its derivatives in s
inline void ph_cumhaz(int dist, double t, double s, double& A, double& As, double& Ass) {
  if (!(t > 0.0)) { A = As = Ass = 0.0; return; }
  if (dist == PD_EXPONENTIAL) { A = t; As = Ass = 0.0; return; }
  if (dist == PD_WEIBULL) {
    const double k = std::exp(s), ky = k * std::log(t);
    A = std::exp(ky); As = ky * A; Ass = ky * A * (1.0 + ky);
    return;
  }
  double f, f1, f2; phi_fun(s * t, f, f1, f2);                  // Gompertz: H0 = t phi(a t)
  A = t * f; As = t * t * f1; Ass = t * t * t * f2;
}

// log baseline hazard log h0(t) of a PH family with its derivatives in s
inline void ph_loghaz(int dist, double t, double s, double& L, double& Ls, double& Lss) {
  if (dist == PD_EXPONENTIAL) { L = Ls = Lss = 0.0; return; }
  if (dist == PD_WEIBULL) {
    const double k = std::exp(s), y = std::log(t);
    L = s + (k - 1.0) * y; Ls = 1.0 + k * y; Lss = k * y;
    return;
  }
  L = s * t; Ls = t; Lss = 0.0;
}

}  // namespace par_detail

// contribution of one subject; entry <= 0 means no truncation
inline Loc par_local(const ParSpec& sp, double entry, double t, int d, double u, double s) {
  using namespace par_detail;
  Loc o;
  if (sp.aft == 0) {
    const double E = std::exp(u);
    double A, As, Ass, Ae = 0, Aes = 0, Aess = 0;
    ph_cumhaz(sp.dist, t, s, A, As, Ass);
    if (entry > 0.0) ph_cumhaz(sp.dist, entry, s, Ae, Aes, Aess);
    const double dA = A - Ae, dAs = As - Aes, dAss = Ass - Aess;
    o.l = -E * dA; o.lu = -E * dA; o.luu = -E * dA; o.lus = -E * dAs; o.ls = -E * dAs; o.lss = -E * dAss;
    if (d) {
      double L, Ls, Lss; ph_loghaz(sp.dist, t, s, L, Ls, Lss);
      o.l += u + L; o.lu += 1.0; o.ls += Ls; o.lss += Lss;
    }
    if (!sp.has_shape()) { o.ls = o.lus = o.lss = 0.0; }
    return o;
  }
  const double sigma = sp.has_shape() && sp.dist != PD_EXPONENTIAL ? std::exp(s) : 1.0;
  const double y = std::log(t), w = (y - u) / sigma;
  double g, g1, g2, G, G1, G2;
  aft_error(sp.dist, w, g, g1, g2, G, G1, G2);
  if (d) { aft_term(g, g1, g2, w, sigma, o, 1.0); o.l -= std::log(sigma) + y; o.ls -= 1.0; }
  else aft_term(G, G1, G2, w, sigma, o, 1.0);
  if (entry > 0.0) {
    const double we = (std::log(entry) - u) / sigma;
    double ge, ge1, ge2, Ge, Ge1, Ge2;
    aft_error(sp.dist, we, ge, ge1, ge2, Ge, Ge1, Ge2);
    aft_term(Ge, Ge1, Ge2, we, sigma, o, -1.0);
  }
  if (sp.dist == PD_EXPONENTIAL) { o.ls = o.lus = o.lss = 0.0; }
  return o;
}

// the data of one fit; X row-major n x p, entry may be null
struct ParData {
  int n = 0;
  const double* entry = nullptr;
  const double* time = nullptr;
  const int* status = nullptr;
  const double* X = nullptr;
};

inline double par_u(const ParSpec& sp, const double* xrow, const std::vector<double>& th) {
  double u = th[0];
  for (int j = 0; j < sp.p; ++j) u += xrow[j] * th[1 + j];
  return u;
}

struct ParModel {
  ParSpec sp;
  ParData dat;
  void eval(const std::vector<double>& th, Eval& ev, bool hess) const {
    const int P = sp.npar(), p = sp.p;
    const int is = sp.has_shape() ? P - 1 : -1;
    const double s = is >= 0 ? th[is] : 0.0;
    ev.ll = 0.0; ev.g.assign(P, 0.0);
    if (hess) { ev.H = Mat(P); }
    std::vector<double> x(p + 1);
    for (int i = 0; i < dat.n; ++i) {
      const double* xr = dat.X + static_cast<std::size_t>(i) * p;
      x[0] = 1.0; for (int j = 0; j < p; ++j) x[j + 1] = xr[j];
      const double u = par_u(sp, xr, th);
      const Loc o = par_local(sp, dat.entry ? dat.entry[i] : 0.0, dat.time[i], dat.status[i], u, s);
      ev.ll += o.l;
      for (int a = 0; a <= p; ++a) ev.g[a] += o.lu * x[a];
      if (is >= 0) ev.g[is] += o.ls;
      if (hess) {
        for (int a = 0; a <= p; ++a) {
          for (int b = 0; b <= a; ++b) ev.H(a, b) -= o.luu * x[a] * x[b];
          if (is >= 0) ev.H(is, a) -= o.lus * x[a];
        }
        if (is >= 0) ev.H(is, is) -= o.lss;
      }
    }
    if (!std::isfinite(ev.ll)) ev.ll = -std::numeric_limits<double>::infinity();
    if (hess) for (int a = 0; a < P; ++a) for (int b = a + 1; b < P; ++b) ev.H(a, b) = ev.H(b, a);
  }
};

// Starting values from the exponential model, whose maximum likelihood estimate is closed form:
// rate = events / total time at risk. The shape starts at the value that reduces each family to
// (or near) the exponential.
inline std::vector<double> par_start(const ParSpec& sp, const ParData& dat) {
  double D = 0, T = 0;
  for (int i = 0; i < dat.n; ++i) { D += dat.status[i] ? 1.0 : 0.0; T += dat.time[i] - (dat.entry ? dat.entry[i] : 0.0); }
  const double rate = std::max(D, 0.5) / std::max(T, 1e-300);
  std::vector<double> th(sp.npar(), 0.0);
  if (sp.aft == 0) th[0] = std::log(rate);
  else if (sp.dist == PD_EXPONENTIAL || sp.dist == PD_WEIBULL) th[0] = -std::log(rate);
  else th[0] = -std::log(rate) + std::log(std::log(2.0));     // median of the exponential
  return th;
}

// Newton-Raphson with step halving; when the information is not positive definite (far from the
// maximum), a Levenberg-Marquardt ridge is added to it (More 1978, Lecture Notes in Math 630:105-116).
inline FitResult par_newton(const ParModel& m, std::vector<double> th, int max_iter, double tol) {
  FitResult r;
  const std::size_t P = th.size();
  Eval ev, evt;
  m.eval(th, ev, true);
  r.ll_init = ev.ll;
  if (!std::isfinite(ev.ll)) { r.flags |= FIT_NONFINITE; r.beta = th; r.ll = ev.ll; return r; }
  std::vector<double> step(P), trial(P);
  for (int it = 1; it <= max_iter; ++it) {
    r.iterations = it;
    double lambda = 0.0;
    bool ok = false;
    for (int k = 0; k < 30 && !ok; ++k) {
      Mat L = ev.H;
      if (lambda > 0) for (std::size_t j = 0; j < P; ++j) L(j, j) += lambda * std::fmax(std::fabs(ev.H(j, j)), 1.0);
      if (cholesky(L)) { step = ev.g; chol_solve(L, step); ok = true; }
      else lambda = lambda == 0.0 ? 1e-6 : lambda * 10.0;
    }
    if (!ok) { r.flags |= FIT_INFO_NOT_PD; break; }
    double scale = 1.0; bool accepted = false;
    for (int h = 0; h < 50; ++h) {
      for (std::size_t j = 0; j < P; ++j) trial[j] = th[j] + scale * step[j];
      m.eval(trial, evt, true);
      if (std::isfinite(evt.ll) && evt.ll >= ev.ll - 1e-13 * (std::fabs(ev.ll) + 1.0)) { accepted = true; break; }
      scale *= 0.5;
    }
    if (!accepted) { r.flags |= FIT_NONFINITE; break; }
    const double ll_old = ev.ll;
    r.last_step.assign(P, 0.0);
    for (std::size_t j = 0; j < P; ++j) r.last_step[j] = scale * step[j];
    th = trial;
    std::swap(ev, evt);
    double gmax = 0.0;
    for (double gj : ev.g) gmax = std::fmax(gmax, std::fabs(gj));
    if (lambda == 0.0 && std::fabs(ev.ll - ll_old) <= tol * (std::fabs(ev.ll) + tol) && gmax < 1e-3) { r.converged = true; break; }
  }
  if (!r.converged && r.flags == FIT_OK) r.flags |= FIT_MAX_ITER;
  r.beta = th; r.ll = ev.ll;
  for (double gj : ev.g) r.grad_maxabs = std::fmax(r.grad_maxabs, std::fabs(gj));
  Mat L = ev.H;
  if (cholesky(L)) r.var = chol_inverse(L); else r.flags |= FIT_INFO_NOT_PD;
  return r;
}

// ------------------------------------------------------------------------------------------------
// quantities of a fitted model at covariate row x
// ------------------------------------------------------------------------------------------------
struct ParQ { double H, logh; };           // cumulative hazard and log hazard

inline ParQ par_hazard(const ParSpec& sp, const std::vector<double>& th, const double* xrow, double t) {
  using namespace par_detail;
  const double u = par_u(sp, xrow, th);
  const double s = sp.has_shape() ? th[sp.npar() - 1] : 0.0;
  ParQ q;
  if (!(t > 0.0)) {                       // at time zero: H = 0, and the hazard is its limit as t -> 0
    q.H = 0.0;
    const double k = sp.aft == 0 ? std::exp(s) : (sp.dist == PD_EXPONENTIAL ? 1.0 : std::exp(-s));   // power of t near zero, minus one
    if (sp.aft == 0 && sp.dist != PD_WEIBULL) q.logh = u;                       // exponential, Gompertz: h(0) = exp(eta)
    else if (sp.aft == 1 && sp.dist == PD_LOGNORMAL) q.logh = -INFINITY;        // lognormal hazard starts at zero
    else if (k < 1.0) q.logh = INFINITY;
    else if (k > 1.0) q.logh = -INFINITY;
    else q.logh = sp.aft == 0 ? u + s : -u;                                     // shape exactly one: exponential
    return q;
  }
  if (sp.aft == 0) {
    double A, As, Ass, L, Ls, Lss;
    ph_cumhaz(sp.dist, t, s, A, As, Ass); ph_loghaz(sp.dist, t, s, L, Ls, Lss);
    q.H = std::exp(u) * A; q.logh = u + L;
    return q;
  }
  const double sigma = sp.dist == PD_EXPONENTIAL ? 1.0 : std::exp(s), y = std::log(t), w = (y - u) / sigma;
  double g, g1, g2, G, G1, G2;
  aft_error(sp.dist, w, g, g1, g2, G, G1, G2);
  q.H = -G; q.logh = g - G - std::log(sigma) - y;
  return q;
}

// inverse standard normal distribution function: Acklam's rational approximation refined by one
// Halley step (Acklam 2003; the refinement gives full double precision).
inline double qnorm_std(double p) {
  if (!(p > 0.0)) return -INFINITY;
  if (!(p < 1.0)) return INFINITY;
  static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
  static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01};
  static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
  static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00};
  double x;
  if (p < 0.02425) { const double q = std::sqrt(-2 * std::log(p)); x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  else if (p > 1 - 0.02425) { const double q = std::sqrt(-2 * std::log1p(-p)); x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1); }
  else { const double q = p - 0.5, r = q * q; x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1); }
  const double e = 0.5 * std::erfc(-x / std::sqrt(2.0)) - p;
  const double u = e * 2.5066282746310002 * std::exp(x * x / 2);
  return x - u / (1 + x * u / 2);
}

// time at which the survival function equals 1 - prob (Inf if never reached)
inline double par_quantile(const ParSpec& sp, const std::vector<double>& th, const double* xrow, double prob) {
  const double u = par_u(sp, xrow, th);
  const double s = sp.has_shape() ? th[sp.npar() - 1] : 0.0;
  const double Hq = -std::log1p(-prob);
  if (sp.aft == 0) {
    const double h = Hq / std::exp(u);
    if (sp.dist == PD_EXPONENTIAL) return h;
    if (sp.dist == PD_WEIBULL) return std::pow(h, std::exp(-s));
    if (std::fabs(s) < 1e-12) return h;
    const double a = 1.0 + s * h;
    return a > 0.0 ? std::log(a) / s : INFINITY;
  }
  const double sigma = sp.dist == PD_EXPONENTIAL ? 1.0 : std::exp(s);
  double w;
  if (sp.dist == PD_EXPONENTIAL || sp.dist == PD_WEIBULL) w = std::log(Hq);
  else if (sp.dist == PD_LOGNORMAL) w = qnorm_std(prob);
  else w = std::log(prob / (1.0 - prob));
  return std::exp(u + sigma * w);
}

// ------------------------------------------------------------------------------------------------
// numerical integration: Gauss-Legendre with 20 nodes on [0, 1] (nodes and weights to 16 digits)
// ------------------------------------------------------------------------------------------------
namespace par_detail {
inline const double* gl20_nodes() {
  static const double x[10] = {0.0765265211334973, 0.2277858511416451, 0.3737060887154195, 0.5108670019508271, 0.6360536807265150,
                               0.7463319064601508, 0.8391169718222188, 0.9122344282513259, 0.9639719272779138, 0.9931285991850949};
  return x;
}
inline const double* gl20_weights() {
  static const double w[10] = {0.1527533871307258, 0.1491729864726037, 0.1420961093183820, 0.1316886384491766, 0.1181945319615184,
                               0.1019301198172404, 0.0832767415767048, 0.0626720483341091, 0.0406014298003869, 0.0176140071391521};
  return w;
}
// integral of f over [a, b] by 20-point Gauss-Legendre
template <class F> double gl20(F&& f, double a, double b) {
  const double* x = gl20_nodes(); const double* w = gl20_weights();
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  double s = 0.0;
  for (int i = 0; i < 10; ++i) s += w[i] * (f(c - h * x[i]) + f(c + h * x[i]));
  return s * h;
}
// integral over [0, b] of a function with an integrable singularity u^(k-1), k > 0, at zero (a
// Weibull hazard with shape below one): geometrically graded pieces [b 4^-(j+1), b 4^-j], j < 26,
// and the substitution u = c v^3 on the last piece [0, c]
template <class F> double gl20_origin(F&& f, double b) {
  double s = 0.0, hi = b;
  for (int j = 0; j < 26; ++j) { const double lo = hi * 0.25; s += gl20(f, lo, hi); hi = lo; }
  const double c = hi;
  s += gl20([&](double v) { const double u = c * v * v * v; return f(u) * 3.0 * c * v * v; }, 0.0, 1.0);
  return s;
}
}  // namespace par_detail

// Cumulative incidence F_k(t) = int_0^t h_k(u) S(u) du, S(u) = exp(-sum_j H_j(u)), of every cause
// at increasing times t_1 < ... < t_m for one covariate row (Prentice et al. 1978; Kalbfleisch &
// Prentice 2002, Sec. 8.2). Each interval between consecutive times is cut into pieces of length at
// most tmax / 50 and integrated by 20-point Gauss-Legendre, all causes at the same nodes. When a
// hazard is unbounded at zero (Weibull shape below one) the first piece is graded geometrically
// towards zero. out: m x K row-major; surv: m.
inline bool par_singular_at_zero(const ParSpec& sp, const std::vector<double>& th) {
  if (!sp.has_shape()) return false;
  const double s = th[sp.npar() - 1];
  if (sp.aft == 0) return sp.dist == PD_WEIBULL && s < 0.0;
  return (sp.dist == PD_WEIBULL || sp.dist == PD_LOGLOGISTIC) && s > 0.0;
}
inline void par_cif(int K, const ParSpec* sp, const std::vector<double>* th, const double* const* xrow,
                    int m, const double* times, double tmax, double* out, double* surv) {
  std::vector<double> acc(K, 0.0), hz(K);
  bool singular = false;
  for (int k = 0; k < K; ++k) singular = singular || par_singular_at_zero(sp[k], th[k]);
  const double* gx = par_detail::gl20_nodes(); const double* gw = par_detail::gl20_weights();
  // add w * h_k(u) S(u) for every cause
  auto node = [&](double u, double w) {
    if (!(u > 0.0) || w == 0.0) return;
    double H = 0.0;
    for (int k = 0; k < K; ++k) { const ParQ q = par_hazard(sp[k], th[k], xrow[k], u); H += q.H; hz[k] = q.logh; }
    const double S = std::exp(-H);
    if (S == 0.0) return;
    for (int k = 0; k < K; ++k) acc[k] += w * std::exp(hz[k]) * S;
  };
  auto piece = [&](double a, double b) {
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    for (int i = 0; i < 10; ++i) { node(c - h * gx[i], h * gw[i]); node(c + h * gx[i], h * gw[i]); }
  };
  auto origin = [&](double b) {                       // [0, b] graded: [b 4^-(j+1), b 4^-j], then u = c v^3
    double hi = b;
    for (int j = 0; j < 26; ++j) { const double lo = hi * 0.25; piece(lo, hi); hi = lo; }
    const double c = hi;
    for (int i = 0; i < 10; ++i) for (int sgn = -1; sgn <= 1; sgn += 2) {
      const double v = 0.5 + 0.5 * sgn * gx[i];
      node(c * v * v * v, 0.5 * gw[i] * 3.0 * c * v * v);
    }
  };
  auto Ssum = [&](double u) { double H = 0; for (int j = 0; j < K; ++j) H += par_hazard(sp[j], th[j], xrow[j], u).H; return std::exp(-H); };
  const double plen = tmax > 0 ? tmax / 50.0 : 1.0;
  double prev = 0.0;
  for (int i = 0; i < m; ++i) {
    const double t = times[i];
    if (t > prev) {
      const int np = std::max(1, static_cast<int>(std::ceil((t - prev) / plen)));
      const double len = (t - prev) / np;
      for (int q = 0; q < np; ++q) {
        const double a = prev + q * len, b = a + len;
        if (a == 0.0 && singular) origin(b); else piece(a, b);
      }
      prev = t;
    }
    for (int k = 0; k < K; ++k) out[static_cast<std::size_t>(i) * K + k] = acc[k];
    if (surv) surv[i] = Ssum(t);
  }
}

// restricted mean survival time int_0^tau S(u) du of one model at row x
inline double par_rmst(const ParSpec& sp, const std::vector<double>& th, const double* xrow, double tau) {
  if (!(tau > 0)) return 0.0;
  const int np = 50;
  double s = 0.0;
  auto f = [&](double u) { return std::exp(-par_hazard(sp, th, xrow, u).H); };
  for (int q = 0; q < np; ++q) s += par_detail::gl20(f, tau * q / np, tau * (q + 1) / np);
  return s;
}

// Standard errors by the delta method with a central-difference gradient of the predicted quantity
// in the parameters: se^2 = g' V g. q(theta) returns the quantity; h_j = 1e-5 max(1, |theta_j|).
template <class Q> double delta_se(Q&& q, std::vector<double> th, const std::vector<double>& V) {
  const std::size_t P = th.size();
  std::vector<double> g(P);
  for (std::size_t j = 0; j < P; ++j) {
    const double h = 1e-5 * std::fmax(1.0, std::fabs(th[j])), keep = th[j];
    th[j] = keep + h; const double up = q(th);
    th[j] = keep - h; const double dn = q(th);
    th[j] = keep;
    g[j] = (up - dn) / (2.0 * h);
  }
  double v = 0.0;
  for (std::size_t a = 0; a < P; ++a) for (std::size_t b = 0; b < P; ++b) v += g[a] * V[a * P + b] * g[b];
  return v >= 0.0 ? std::sqrt(v) : std::numeric_limits<double>::quiet_NaN();
}

}  // namespace survcore
